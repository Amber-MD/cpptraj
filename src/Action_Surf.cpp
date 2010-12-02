// Distance
#include <cstdio>
#include <cstdlib>
#include "Action_Surf.h"
#include "vectormath.h" // For FOURPI

// CONSTRUCTOR
Surf::Surf() {
  //fprintf(stderr,"Surf Con\n");
  surf=NULL;
  distances=NULL;
  soluteAtoms=0;
} 

// DESTRUCTOR
Surf::~Surf() { 
  if (distances!=NULL) free(distances);
}

/*
 * Surf::init()
 * Expected call: surf <name> <mask1> [out filename]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Surf::init() {
  char *mask1;
  char *surfFile;

  // Get keywords
  surfFile = A->getKeyString("out",NULL);

  // Get Masks
  mask1 = A->getNextMask();
  Mask1.SetMaskString(mask1);

  // Dataset to store surface area 
  surf = DSL->Add(DOUBLE, A->getNextString(),"SA");
  if (surf==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(surfFile,surf);

  //dist->Info();
  fprintf(stdout,"    SURF: %s\n",Mask1.maskString);

  return 0;
}

/*
 * Surf::setup()
 * Set surf up for this parmtop. Get masks etc.
 * P is set in Action::Setup
 */
int Surf::setup() {
  int matrixSize, i;

  if ( Mask1.SetupMask(P,debug) ) return 1;
  if (Mask1.None()) {
    fprintf(stdout,"    Error: Surf::setup: Mask contains 0 atoms.\n");
    return 1;
  }

  // Setup surface area calc for this parm
  P->SetSurfaceInfo();

  // If solvent is present figure out how many solute atoms
  if (P->firstSolvMol > 0) {
    soluteAtoms = 0;
    i = 0;
    while (i < P->firstSolvMol) soluteAtoms += P->atomsPerMol[i++];
  } else {
    soluteAtoms = P->natom;
  }
  fprintf(stdout,"    SURF: Setting up %i solute atoms for surface calc.\n",soluteAtoms);

  // Allocate distance array
  // Use half square matrix minus the diagonal
  matrixSize = ( (soluteAtoms * soluteAtoms) - soluteAtoms ) / 2;
  distances = (double*) realloc(distances, matrixSize * sizeof(double));

  // Reserve space for neighbor list
  NeighborList.reserve(soluteAtoms);

  return 0;  
}

/* 
 * CalcIndex()
 * Given a row and column value for a half square matrix minus diagonal 
 * calculate the position in a 1D array.
 * NOTE: Is all the extra math in index calc worth it? Just have simple
 * square matrix?
 */
int CalcIndex(int iIn, int jIn, int natom) {
  int i, j, i1;

  if (iIn > jIn) {
    j = iIn;
    i = jIn;
  } else {
    i = iIn;
    j = jIn;
  }

  i1 = i + 1;
  return ( ( (natom * i) - ((i1 * i) / 2) ) + j - i1 );
}

/*
 * Surf::CalcLCPO()
 * Return the surface area for the given atom.
 */

/*
 * Surf::action()
 * Calculate surface area.
 * LCPO method from:
 *   J. Weiser, P.S. Shenkin, and W.C. Still,
 *   "Approximate atomic surfaces from linear combinations of pairwise
 *   overlaps (LCPO)", J. Comp. Chem. 20:217 (1999).
 * Adapted from gbsa=1 method in SANDER, egb.f and gbsa.h
 */
int Surf::action() {
  double SA; // DEBUG
  double Si,vdwi,vdwj,vdwk,vdwi2,vdwj2;
  double dij,djk,tmpaij,tmpajk,vdw2dif,ai,aij,ajk; // NOTE: vdw2dif necessary?
  double sumaij, sumajk, sumajk_2, sumaijajk;
  int atomi, atomj, distIndex;
  std::vector<int> ineighbor;
  std::vector<int>::iterator jt;
  std::vector<int>::iterator kt;

  SA = 0.0;

  // Step 1: Calculate all distances 
  distIndex = 0;
  for (atomi = 0; atomi < soluteAtoms - 1; atomi++) {
    for (atomj = atomi + 1; atomj < soluteAtoms; atomj++) {
      distances[distIndex++] = F->DIST(atomi, atomj);
      // DEBUG
      //fprintf(stdout,"SURF_DIST:  %i %i %i %lf\n",atomi,atomj,distIndex-1,distances[distIndex-1]);
    }
  } 
      
  // Step 2: Set up neighbor list for each atom
  NeighborList.clear(); 
  for (atomi = 0; atomi < soluteAtoms; atomi++) {
    // Set up neighbor list for atom i
    ineighbor.clear();
    for (atomj = 0; atomj < soluteAtoms; atomj++) {
      if (atomi==atomj) continue;
      distIndex = CalcIndex(atomi, atomj, soluteAtoms);
      // DEBUG
      //fprintf(stdout,"SURF_NEIG:  %i %i %i %lf\n",atomi,atomj,distIndex,distances[distIndex]);
      // Count atoms as neighbors if their VDW radii touch
      if ( (P->SurfaceInfo[atomi].vdwradii + P->SurfaceInfo[atomj].vdwradii) >
           distances[distIndex] ) {
        // Consider all atoms for icosa, only non-H's for LCPO
        if ( P->SurfaceInfo[atomi].vdwradii > 2.5 &&
             P->SurfaceInfo[atomj].vdwradii > 2.5 ) {
          ineighbor.push_back(atomj);
        }
      }
    }
    NeighborList.push_back(ineighbor);
  }  

  // Step 3: Calculate LCPO surface area
  // Ai = P1*S1 + P2*Sum(Aij) + P3*Sum(Ajk) + P4*Sum(Aij * Sum(Ajk))
  for (atomi = 0; atomi < soluteAtoms; atomi++) {
    sumaij = 0.0;
    sumajk = 0.0;
    sumaijajk = 0.0;

    // DEBUG - print neighbor list
    //fprintf(stdout,"SURF: Neighbors for atom %i:",atomi);
    //ineighbor = NeighborList[atomi];
    //for (jt = ineighbor.begin(); jt != ineighbor.end(); jt++) {
    //  fprintf(stdout," %i",*jt);
    //}
    //fprintf(stdout,"\n");

    // Calculate surface area of atom i
    vdwi = P->SurfaceInfo[atomi].vdwradii;
    vdwi2 = vdwi * vdwi;
    Si = vdwi2 * FOURPI;
    
    // Loop over all neighbors of atomi (j)
    // NOTE: Factor through the 2 in aij and ajk?
    ineighbor = NeighborList[atomi];
    for (jt = ineighbor.begin(); jt != ineighbor.end(); jt++) {
      //printf("i,j %i %i\n",atomi + 1,(*jt)+1);
      vdwj = P->SurfaceInfo[*jt].vdwradii;
      vdwj2 = vdwj * vdwj;
      distIndex = CalcIndex(atomi, *jt, soluteAtoms);
      dij = distances[distIndex];
      tmpaij = vdwi - (dij * 0.5) - ( (vdwi2 - vdwj2)/(2.0 * dij) );
      aij = 2.0 * PI * vdwi * tmpaij;
      sumaij += aij;

      // Find which neighbors of atom i (j and k) are themselves neighbors
      sumajk_2 = 0.0;
      for (kt = ineighbor.begin(); kt != ineighbor.end(); kt++) {
        if ( (*kt) == (*jt) ) continue;
        //printf("i,j,k %i %i %i\n",atomi + 1,(*jt)+1,(*kt)+1);
        vdwk = P->SurfaceInfo[*kt].vdwradii;
        distIndex = CalcIndex(*jt, *kt, soluteAtoms);
        djk = distances[distIndex];
        //printf("%4s%6i%6i%12.8lf\n","DJK ",(*jt)+1,(*kt)+1,djk);
        //printf("%6s%6.2lf%6.2lf\n","AVD ",vdwj,vdwk);
        if ( (vdwj + vdwk) > djk ) {
          vdw2dif = vdwj2 - (vdwk * vdwk);
          tmpajk = (2.0*vdwj) - djk - (vdw2dif / djk);
          ajk = PI*vdwj*tmpajk;
          //tmpajk = vdwj - (djk *0.5) - ( (vdwj2 - (vdwk * vdwk))/(2.0 * djk) );
          //ajk = 2.0 * PI * vdwi * tmpajk;
          //printf("%4s%6i%6i%12.8lf%12.8lf%12.8lf\n","AJK ",(*jt)+1,(*kt)+1,ajk,vdw2dif,tmpajk);
          sumajk += ajk;
          sumajk_2 += ajk;
        }
      } // END loop over neighbor-neighbor pairs of atom i (k)

      sumaijajk += (aij * sumajk_2);

      // DEBUG
      //printf("%4s%20.8lf %20.8lf %20.8lf\n","AJK ",aij,sumajk,sumaijajk);

    } // END Loop over neighbors of atom i (j)

    ai = (P->SurfaceInfo[atomi].P1 * Si) +
         (P->SurfaceInfo[atomi].P2 * sumaij) +
         (P->SurfaceInfo[atomi].P3 * sumajk) +
         (P->SurfaceInfo[atomi].P4 * sumaijajk);
    SA += ai;
  }

  surf->Add(currentFrame, &SA);

  return 0;
} 



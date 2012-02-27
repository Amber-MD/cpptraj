// Surf 
#include "Action_Surf.h"
#include "Constants.h" // For FOURPI, TWOPI
#include "CpptrajStdio.h"
#ifdef _OPENMP
#  include "omp.h"
#endif
// CONSTRUCTOR
Surf::Surf() {
  //fprintf(stderr,"Surf Con\n");
  surf=NULL;
} 

// Surf::init()
/** Expected call: surf <name> <mask1> [out filename]
  */
int Surf::init() {
  char *mask1;
  char *surfFile;

  // Get keywords
  surfFile = actionArgs.getKeyString("out",NULL);

  // Get Masks
  mask1 = actionArgs.getNextMask();
  Mask1.SetMaskString(mask1);

  // Dataset to store surface area 
  surf = DSL->Add(DOUBLE, actionArgs.getNextString(),"SA");
  if (surf==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(surfFile,surf);

  mprintf("    SURF: Calculating surface area for atoms in mask [%s]\n",Mask1.MaskString());

  return 0;
}

// Surf::setup()
/** Set LCPO surface area calc parameters for this parmtop if not already set. 
  * Get the mask, and check that the atoms in mask belong to solute. 
  */
int Surf::setup() {
  int soluteAtoms;

  if (currentParm->SetupIntegerMask( Mask1, activeReference)) return 1;
  if (Mask1.None()) {
    mprintf("    Error: Surf::setup: Mask contains 0 atoms.\n");
    return 1;
  }

  // Setup surface area calc for this parm
  soluteAtoms = currentParm->SetSurfaceInfo();
  if (soluteAtoms < 0) {
    mprinterr("Error: Surf: Could not set up surface parameters for %s.\n",currentParm->parmName);
    return 1;
  }
  mprintf("      SURF: Set up parameters for %i solute atoms.\n",soluteAtoms);
  mprintf("            LCPO surface area will be calculated for %i atoms.\n",Mask1.Nselected);

  // Check that each atom in Mask1 is part of solute
  // Create a separate mask for building the atom neighbor list
  // that only includes atoms with vdw radius > 2.5
  // Consider all atoms for icosa, only non-H's for LCPO
  atomi_neighborMask.ResetMask();
  atomi_noNeighborMask.ResetMask();
  atomj_neighborMask.ResetMask();
  for (int i=0; i < Mask1.Nselected; i++) {
    int atomi = Mask1.Selected[i];
    if (currentParm->SurfaceInfo[atomi].vdwradii > 2.5)
      atomi_neighborMask.AddAtom(atomi);
    else
      atomi_noNeighborMask.AddAtom(atomi);
    if (atomi >= soluteAtoms) {
      mprintf("Error: Surf::setup(): Atom %i in mask %s does not belong to solute.\n",
              atomi+1, Mask1.MaskString());
      return 1;
    }
  }
  // From all atoms, create a second mask for building atom neighbor list
  // that only includes atoms with vdw radius > 2.5
  for (int atomj=0; atomj < soluteAtoms; atomj++) {
    if (currentParm->SurfaceInfo[atomj].vdwradii > 2.5)
      atomj_neighborMask.AddAtom(atomj);
  }
 
  return 0;  
}

// Surf::action()
/** Calculate surface area. */
int Surf::action() {
  double SA; 
  int atomi, atomj;
  std::vector<int> ineighbor;
  std::vector<double> Distances_i_j;

  // Set up neighbor list for each atom in mask and calc its surface 
  // area contribution. Sum these up to get the total surface area.
  SA = 0.0;
#ifdef _OPENMP
#pragma omp parallel private(atomi,atomj,ineighbor,Distances_i_j) reduction(+: SA)
{
#pragma omp for 
#endif      
  for (int maskIndex = 0; maskIndex < atomi_neighborMask.Nselected; maskIndex++) {
    atomi = atomi_neighborMask.Selected[maskIndex];
    // Vdw of atom i
    double vdwi = currentParm->SurfaceInfo[atomi].vdwradii;
    // Set up neighbor list for atom i
    ineighbor.clear();
    Distances_i_j.clear();
    for (int ajidx = 0; ajidx < atomj_neighborMask.Nselected; ajidx++) {
      atomj = atomj_neighborMask.Selected[ajidx];
      if (atomi!=atomj) {
        double dij = currentFrame->DIST(atomi, atomj);
        // Count atoms as neighbors if their VDW radii touch
        if ( (vdwi + currentParm->SurfaceInfo[atomj].vdwradii) > dij ) {
          ineighbor.push_back(atomj);
          Distances_i_j.push_back(dij);
        }
        //mprintf("SURF_NEIG:  %i %i %lf\n",atomi,atomj,dij);
      }
    }
    // Calculate surface area
    // -------------------------------------------------------------------------
    // Calculate LCPO surface area Ai for atomi:
    // Ai = P1*S1 + P2*Sum(Aij) + P3*Sum(Ajk) + P4*Sum(Aij * Sum(Ajk))
    double sumaij = 0.0;
    double sumajk = 0.0;
    double sumaijajk = 0.0;

    // DEBUG - print neighbor list
    //mprintf("SURF: Neighbors for atom %i:",atomi);
    //for (std::vector<int>::iterator jt = ineighbor.begin(); jt != ineighbor.end(); jt++) {
    //  mprintf(" %i",*jt);
    //}
    //mprintf("\n");

    // Calculate surface area of atom i
    double vdwi2 = vdwi * vdwi;
    double Si = vdwi2 * FOURPI;

    // Loop over all neighbors of atomi (j)
    // NOTE: Factor through the 2 in aij and ajk?
    std::vector<double>::iterator Dij = Distances_i_j.begin();
    for (std::vector<int>::iterator jt = ineighbor.begin(); 
                                    jt != ineighbor.end(); 
                                    jt++) 
    {
      //printf("i,j %i %i\n",atomi + 1,(*jt)+1);
      double vdwj = currentParm->SurfaceInfo[*jt].vdwradii;
      double vdwj2 = vdwj * vdwj;
      double dij = *Dij;
      double tmpaij = vdwi - (dij * 0.5) - ( (vdwi2 - vdwj2)/(2.0 * dij) );
      double aij = TWOPI * vdwi * tmpaij;
      sumaij += aij;

      // Find which neighbors of atom i (j and k) are themselves neighbors
      double sumajk_2 = 0.0;
      for (std::vector<int>::iterator kt = ineighbor.begin(); 
                                      kt != ineighbor.end(); 
                                      kt++) 
      {
        if ( (*kt) == (*jt) ) continue;
        //printf("i,j,k %i %i %i\n",atomi + 1,(*jt)+1,(*kt)+1);
        double vdwk = currentParm->SurfaceInfo[*kt].vdwradii;
        double djk = currentFrame->DIST(*jt, *kt);
        //printf("%4s%6i%6i%12.8lf\n","DJK ",(*jt)+1,(*kt)+1,djk);
        //printf("%6s%6.2lf%6.2lf\n","AVD ",vdwj,vdwk);
        if ( (vdwj + vdwk) > djk ) {
          double vdw2dif = vdwj2 - (vdwk * vdwk);
          double tmpajk = (2.0*vdwj) - djk - (vdw2dif / djk);
          double ajk = PI*vdwj*tmpajk;
          //tmpajk = vdwj - (djk *0.5) - ( (vdwj2 - (vdwk * vdwk))/(2.0 * djk) );
          //ajk = 2.0 * PI * vdwi * tmpajk;
          //printf("%4s%6i%6i%12.8lf%12.8lf%12.8lf\n","AJK ",(*jt)+1,(*kt)+1,ajk,vdw2dif,tmpajk);
          sumajk += ajk;
          sumajk_2 += ajk;
        }
      } // END loop over neighbor-neighbor pairs of atom i (kt)

      sumaijajk += (aij * sumajk_2);

      // DEBUG
      //printf("%4s%20.8lf %20.8lf %20.8lf\n","AJK ",aij,sumajk,sumaijajk);
      Dij++;

    } // END Loop over neighbors of atom i (jt)
    SA += ( (currentParm->SurfaceInfo[atomi].P1 * Si       ) +
            (currentParm->SurfaceInfo[atomi].P2 * sumaij   ) +
            (currentParm->SurfaceInfo[atomi].P3 * sumajk   ) +
            (currentParm->SurfaceInfo[atomi].P4 * sumaijajk) 
          );
    // -------------------------------------------------------------------------
  } // END Loop over atoms in mask (atomi) 
#ifdef _OPENMP
} // END pragma omp parallel
#endif
  // Second loop over atoms with no neighbors (vdw <= 2.5)
  for (int maskIndex=0; maskIndex < atomi_noNeighborMask.Nselected; maskIndex++) {
    int atomi = atomi_noNeighborMask.Selected[maskIndex];
    // Vdw of atom i
    double vdwi = currentParm->SurfaceInfo[atomi].vdwradii;
    // Calculate surface area of atom i
    double vdwi2 = vdwi * vdwi;
    double Si = vdwi2 * FOURPI; 
    SA += (currentParm->SurfaceInfo[atomi].P1 * Si);
  }

  surf->Add(frameNum, &SA);

  return 0;
} 


// Action_Surf 
#include "Action_Surf.h"
#include "Constants.h" // For FOURPI, TWOPI
#include "CpptrajStdio.h"
#ifdef _OPENMP
#  include "omp.h"
#endif
// CONSTRUCTOR
Action_Surf::Action_Surf() {
  //fprintf(stderr,"Surf Con\n");
  surf=NULL;
} 

// Action_Surf::init()
/** Expected call: surf <name> <mask1> [out filename]
  */
int Action_Surf::init() {
  char *mask1;
  char *surfFile;

  // Get keywords
  surfFile = actionArgs.getKeyString("out",NULL);

  // Get Masks
  mask1 = actionArgs.getNextMask();
  Mask1.SetMaskString(mask1);

  // Dataset to store surface area 
  surf = DSL->Add(DataSet::DOUBLE, actionArgs.getNextString(),"SA");
  if (surf==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(surfFile,surf);

  mprintf("    SURF: Calculating surface area for atoms in mask [%s]\n",Mask1.MaskString());

  return 0;
}

// Action_Surf::setup()
/** Set LCPO surface area calc parameters for this parmtop if not already set. 
  * Get the mask, and check that the atoms in mask belong to solute. 
  */
int Action_Surf::setup() {
  SurfInfo SI;
  int soluteAtoms;

  if (currentParm->SetupIntegerMask( Mask1 )) return 1;
  if (Mask1.None()) {
    mprintf("    Error: Surf::setup: Mask contains 0 atoms.\n");
    return 1;
  }

  // Setup surface area calc for this parm
  soluteAtoms = currentParm->SoluteAtoms();
  if (soluteAtoms <= 0) {
    mprinterr("Error: Surf: No solute atoms in %s.\n",currentParm->c_str());
    return 1;
  }
  mprintf("      SURF: %i solute atoms.\n",soluteAtoms);
  mprintf("            LCPO surface area will be calculated for %i atoms.\n",Mask1.Nselected());

  // Check that each atom in Mask1 is part of solute
  // Create a separate mask for building the atom neighbor list
  // that only includes atoms with vdw radius > 2.5
  // Consider all atoms for icosa, only non-H's for LCPO
  atomi_neighborMask.ResetMask();
  atomi_noNeighborMask.ResetMask();
  atomj_neighborMask.ResetMask();
  SurfaceInfo_neighbor.clear();
  SurfaceInfo_noNeighbor.clear();
  for (AtomMask::const_iterator atomi = Mask1.begin(); atomi!=Mask1.end(); atomi++) {
    SetAtomLCPO( *atomi, (*currentParm)[*atomi], &SI ); 
    if (SI.vdwradii > 2.5) {
      atomi_neighborMask.AddAtom(*atomi);
      SurfaceInfo_neighbor.push_back( SI );
    } else {
      atomi_noNeighborMask.AddAtom(*atomi);
      SurfaceInfo_noNeighbor.push_back( SI );
    }
    if (*atomi >= soluteAtoms) {
      mprintf("Error: Surf::setup(): Atom %i in mask %s does not belong to solute.\n",
              *atomi+1, Mask1.MaskString());
      return 1;
    }
  }
  // From all atoms, create a second mask for building atom neighbor list
  // that only includes atoms with vdw radius > 2.5
  VDW.clear();
  VDW.reserve( soluteAtoms );
  for (int atomj=0; atomj < soluteAtoms; atomj++) {
    SetAtomLCPO( atomj, (*currentParm)[atomj], &SI );
    VDW.push_back( SI.vdwradii );
    if (SI.vdwradii > 2.5)
      atomj_neighborMask.AddAtom(atomj);
  }
 
  return 0;  
}

// Action_Surf::action()
/** Calculate surface area. */
int Action_Surf::action() {
  double SA;
  int atomi, idx;
  AtomMask::const_iterator atomj; 
  std::vector<int> ineighbor;
  std::vector<double> Distances_i_j;
  int max_atomi_neighbormask = atomi_neighborMask.Nselected();

  // Set up neighbor list for each atom in mask and calc its surface 
  // area contribution. Sum these up to get the total surface area.
  SA = 0.0;
#ifdef _OPENMP
#pragma omp parallel private(atomi,idx,atomj,ineighbor,Distances_i_j) reduction(+: SA)
{
#pragma omp for 
#endif
  for (idx = 0; idx < max_atomi_neighbormask; idx++) {
    atomi = atomi_neighborMask[idx];
    // Vdw of atom i
    double vdwi = VDW[atomi];
    // Set up neighbor list for atom i
    ineighbor.clear();
    Distances_i_j.clear();
    for (atomj = atomj_neighborMask.begin(); atomj != atomj_neighborMask.end(); atomj++)
    {
      if (atomi != *atomj) {
        double dij = currentFrame->DIST(atomi, *atomj);
        // Count atoms as neighbors if their VDW radii touch
        if ( (vdwi + VDW[*atomj]) > dij ) {
          ineighbor.push_back(*atomj);
          Distances_i_j.push_back(dij);
        }
        //mprintf("SURF_NEIG:  %i %i %lf\n",atomi,*atomj,dij);
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
      double vdwj = VDW[*jt];
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
        double vdwk = VDW[*kt];
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
    SA += ( (SurfaceInfo_neighbor[idx].P1 * Si       ) +
            (SurfaceInfo_neighbor[idx].P2 * sumaij   ) +
            (SurfaceInfo_neighbor[idx].P3 * sumajk   ) +
            (SurfaceInfo_neighbor[idx].P4 * sumaijajk) 
          );
    // -------------------------------------------------------------------------
  } // END Loop over atoms in mask (atomi) 
#ifdef _OPENMP
} // END pragma omp parallel
#endif
  // Second loop over atoms with no neighbors (vdw <= 2.5)
  int idxj=0;
  for (atomj = atomi_noNeighborMask.begin(); atomj != atomi_noNeighborMask.end(); atomj++)
  {
    // Vdw of atom i
    double vdwi = VDW[*atomj];
    // Calculate surface area of atom i
    double vdwi2 = vdwi * vdwi;
    double Si = vdwi2 * FOURPI; 
    SA += (SurfaceInfo_noNeighbor[idxj++].P1 * Si);
  }

  surf->Add(frameNum, &SA);

  return 0;
} 

// -----------------------------------------------------------------------------
// Action_Surf::AssignLCPO()
/** Assign parameters for LCPO method. All radii are incremented by 1.4 Ang.
  */
void Action_Surf::AssignLCPO(SurfInfo *S, double vdwradii, double P1, double P2,
                      double P3, double P4) 
{
  S->vdwradii = vdwradii + 1.4;
  S->P1 = P1;
  S->P2 = P2;
  S->P3 = P3;
  S->P4 = P4;
}

// WarnLCPO()
/// Called when the number of bonds to the atom of type atype is not usual.
static void WarnLCPO(NameType &atype, int atom, int numBonds) {
  mprintf("Warning: Unusual number of bonds for atom %i (%i), type %s.\n",
          atom, numBonds, *atype);
  mprintf("Using default atom parameters.\n");
}

// Action_Surf::SetAtomLCPO()
/** Set up parameters only used in surface area calcs.
  * LCPO method from:
  *   J. Weiser, P.S. Shenkin, and W.C. Still,
  *   "Approximate atomic surfaces from linear combinations of pairwise
  *   overlaps (LCPO)", J. Comp. Chem. 20:217 (1999).
  * Adapted from gbsa=1 method in SANDER, mdread.f
  * \return the number of solute atoms for which paramters were set. 
  * \return -1 on error.
  */
void Action_Surf::SetAtomLCPO(int i, const Atom &atom, SurfInfo *SIptr) {
  // Get the number of non-H bonded neighbors to this atom
  int numBonds = 0;
  for (Atom::bond_iterator batom = atom.bondbegin(); batom != atom.bondend(); batom++)
    if ( (*currentParm)[ *batom ].Element() != Atom::HYDROGEN )
      ++numBonds;
  NameType atype = atom.Type();

  // Check: Only set parameters for solute atoms?
  // Set vdw radii and LCPO parameters for this atom
  if (atype[0]=='C' && atype[1]=='T') {
    switch ( numBonds ) {
      case 1: AssignLCPO(SIptr, 1.70, 0.77887, -0.28063, -0.0012968, 0.00039328); break;
      case 2: AssignLCPO(SIptr, 1.70, 0.56482, -0.19608, -0.0010219, 0.0002658);  break;
      case 3: AssignLCPO(SIptr, 1.70, 0.23348, -0.072627, -0.00020079, 0.00007967); break;
      case 4: AssignLCPO(SIptr, 1.70, 0.00000, 0.00000, 0.00000, 0.00000); break;
      default: WarnLCPO(atype,i,numBonds);
              AssignLCPO(SIptr, 1.70, 0.77887, -0.28063, -0.0012968, 0.00039328);
    }
  } else if (atype[0]=='C' || atype[0]=='c') {
    switch ( numBonds ) {
      case 2: AssignLCPO(SIptr, 1.70, 0.51245, -0.15966, -0.00019781, 0.00016392); break;
      case 3: AssignLCPO(SIptr, 1.70, 0.070344, -0.019015, -0.000022009, 0.000016875); break;
      default: WarnLCPO(atype,i,numBonds);
              AssignLCPO(SIptr, 1.70, 0.77887, -0.28063, -0.0012968, 0.00039328);
    }
  } else if (atype[0]=='O' && atype[1]==' ') {
    AssignLCPO(SIptr, 1.60, 0.68563, -0.1868, -0.00135573, 0.00023743);
  } else if (atype[0]=='O' && atype[1]=='2') {
    AssignLCPO(SIptr, 1.60, 0.88857, -0.33421, -0.0018683, 0.00049372);
  } else if (atype[0]=='O' || atype[0]=='o') {
    switch (numBonds) {
      case 1: AssignLCPO(SIptr, 1.60, 0.77914, -0.25262, -0.0016056, 0.00035071); break;
      case 2: AssignLCPO(SIptr, 1.60, 0.49392, -0.16038, -0.00015512, 0.00016453); break;
      default: WarnLCPO(atype,i,numBonds);
              AssignLCPO(SIptr, 1.60, 0.77914, -0.25262, -0.0016056, 0.00035071);
    }
  } else if (atype[0]=='N' && atype[1]=='3') {
    switch (numBonds) {
      case 1: AssignLCPO(SIptr, 1.65, 0.078602, -0.29198, -0.0006537, 0.00036247); break;
      case 2: AssignLCPO(SIptr, 1.65, 0.22599, -0.036648, -0.0012297, 0.000080038); break;
      case 3: AssignLCPO(SIptr, 1.65, 0.051481, -0.012603, -0.00032006, 0.000024774); break;
      default: WarnLCPO(atype,i,numBonds);
              AssignLCPO(SIptr, 1.65, 0.078602, -0.29198, -0.0006537, 0.00036247);
    }
  } else if (atype[0]=='N' || atype[0]=='n') {
    switch (numBonds) {
      case 1: AssignLCPO(SIptr, 1.65, 0.73511, -0.22116, -0.00089148, 0.0002523); break;
      case 2: AssignLCPO(SIptr, 1.65, 0.41102, -0.12254, -0.000075448, 0.00011804); break;
      case 3: AssignLCPO(SIptr, 1.65, 0.062577, -0.017874, -0.00008312, 0.000019849); break;
      default: WarnLCPO(atype,i,numBonds);
              AssignLCPO(SIptr, 1.65, 0.078602, -0.29198, -0.0006537, 0.00036247);
    }
  } else if (atype[0]=='S' && atype[1]=='H') {
    AssignLCPO(SIptr, 1.90, 0.7722, -0.26393, 0.0010629, 0.0002179);
  } else if (atype[0]=='S' || atype[0]=='s') {
    AssignLCPO(SIptr, 1.90, 0.54581, -0.19477, -0.0012873, 0.00029247);
  } else if (atype[0]=='P' || atype[1]=='p') {
    switch (numBonds) {
      case 3: AssignLCPO(SIptr, 1.90, 0.3865, -0.18249, -0.0036598, 0.0004264); break;
      case 4: AssignLCPO(SIptr, 1.90, 0.03873, -0.0089339, 0.0000083582, 0.0000030381); break;
      default: WarnLCPO(atype,i,numBonds);
        AssignLCPO(SIptr, 1.90, 0.3865, -0.18249, -0.0036598, 0.0004264);
    }
  } else if (atype[0]=='Z') {
    AssignLCPO(SIptr, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000);
  } else if (atype[0]=='H' || atype[0]=='h') {
    AssignLCPO(SIptr, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000);
  } else if (atype[0]=='M' && atype[1]=='G') {
    //  Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
    //  Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
    //  Mg radius = 1.45A: Aqvist 1992
    //  The following P1-4 values were taken from O.sp3 with two bonded 
    //  neighbors -> O has the smallest van der Waals radius 
    //  compared to all other elements which had been parametrized
    AssignLCPO(SIptr, 1.18, 0.49392, -0.16038, -0.00015512, 0.00016453);
  } else {
    mprintf("Warning: Using carbon SA parms for unknown atom type %i %2\n",i,*atype);
    AssignLCPO(SIptr, 1.70, 0.51245, -0.15966, -0.00019781, 0.00016392);
  }

  // DEBUG
  /*
  for (i=0; i<numSoluteAtoms; i++) {
    fprintf(stdout,"%6i %4s: %6.2lf %lf %lf %lf %lf\n",i+1,types[i],SurfaceInfo[i].vdwradii,
    fprintf(stdout,"%6i%6.2lf%12.8lf%12.8lf%12.8lf%12.8lf\n",i+1,SurfaceInfo[i].vdwradii,
            SurfaceInfo[i].P1,SurfaceInfo[i].P2,SurfaceInfo[i].P3,SurfaceInfo[i].P4);
  }
  */
}


#include <cmath>  // sqrt
#include <cctype> // toupper
#include "Action_Surf.h"
#include "Constants.h" // For FOURPI, TWOPI
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#ifdef _OPENMP
#  include <omp.h>
#endif
// CONSTRUCTOR
Action_Surf::Action_Surf() :
  surf_(0),
  neighborCut_(2.5),
  noNeighborTerm_(0.0)
{} 

// Action_Surf::Help()
void Action_Surf::Help() const {
  mprintf("\t<name> [<mask1>] [out <filename>] [solutemask <mask>]\n"
          "  Calculate LCPO surface area of atoms in <mask1>. If 'solutemask'\n"
          "  is specified, calculate the contribution of <mask1> to selected\n"
          "  solute atoms. If solute is not specified, it is considered to be\n"
          "  any molecule not marked as solvent > 1 atom in size.\n");
}

// Action_Surf::Init()
Action::RetType Action_Surf::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get keywords
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs);

  // Get Masks
  std::string maskexp = actionArgs.GetMaskNext();
  if (!maskexp.empty()) Mask1_.SetMaskString( maskexp );
  maskexp = actionArgs.GetStringKey("solutemask");
  if (!maskexp.empty()) SoluteMask_.SetMaskString( maskexp );

  // Dataset to store surface area 
  surf_ = init.DSL().AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "SA");
  if (surf_==0) return Action::ERR;
  // Add DataSet to DataFileList
  if (outfile != 0) outfile->AddDataSet( surf_ );

  mprintf("    SURF: ");
  if (!Mask1_.MaskStringSet())
    mprintf("Calculating LCPO surface area for all solute atoms.\n");
  else
    mprintf("Calculating LCPO surface area for atoms in mask '%s'\n", Mask1_.MaskString());
  if (!SoluteMask_.MaskStringSet())
    mprintf("\tSolute will be all molecules not marked as solvent with size > 1 atom.\n");
  else
    mprintf("\tSolute will be atoms selected by mask '%s'\n", SoluteMask_.MaskString());
  if (outfile != 0) mprintf("\tOutput to '%s'\n", outfile->DataFilename().full());
  mprintf("#Citation: Weiser, J.; Shenkin, P. S.; Still, W. C.; \"Approximate atomic\n"
          "#          surfaces from linear combinations of pairwise overlaps (LCPO).\"\n"
          "#          J. Comp. Chem. (1999), V.20, pp.217-230.\n");

  return Action::OK;
}

// Action_Surf::Setup()
/** There are two modes. If a mask has been specified, want surface area of
  * only atoms in mask. Otherwise want surface area of all solute atoms. In
  * either case attempt to identify solute atoms (i.e. exclude solvent/ions).
  */
Action::RetType Action_Surf::Setup(ActionSetup& setup) {
  if (Mask1_.MaskStringSet()) {
    if (setup.Top().SetupIntegerMask( Mask1_ )) return Action::ERR;
    if (Mask1_.None()) {
      mprintf("Warning: Mask '%s' corresponds to 0 atoms.\n", Mask1_.MaskString());
      return Action::SKIP;
    }
    Mask1_.MaskInfo();
  }
 
  // Determine solute atoms.
  if (SoluteMask_.MaskStringSet()) {
    if (setup.Top().SetupIntegerMask( SoluteMask_ )) return Action::ERR;
    SoluteMask_.MaskInfo();
    if (SoluteMask_.None()) {
      mprintf("Warning: Solute mask selects no atoms.\n");
      return Action::SKIP;
    }
  } else {
    SoluteMask_.ResetMask();
    SoluteMask_.SetNatoms( setup.Top().Natom() );
    if ( setup.Top().Nmol() > 0) {
      mprintf("\tConsidering only non-solvent molecules with size > 1 as solute.\n");
      for (int at = 0; at != setup.Top().Natom(); at++)
      {
        Molecule const& mol = setup.Top().Mol( setup.Top()[at].MolNum() );
        if (!mol.IsSolvent() && mol.NumAtoms() > 1)
          SoluteMask_.AddSelectedAtom( at );
      }
    } else {
      mprintf("Warning: No molecule info in '%s'. Considering all atoms as solute.\n");
      for (int at = 0; at != setup.Top().Natom(); at++)
        SoluteMask_.AddSelectedAtom( at );
    }
  }

  // If no SA mask specified use all solute atoms.
  if (!Mask1_.MaskStringSet())
    Mask1_ = SoluteMask_;
  mprintf("\t%i solute atoms. Calculating LCPO surface area for %i atoms.\n",
          SoluteMask_.Nselected(), Mask1_.Nselected());
  CharMask cmask1( Mask1_.ConvertToCharMask(), Mask1_.Nselected() );

  // Store vdW radii for all solute atoms with radii greater than neighborCut_.
  // If that atom will have its SA contribution calculated also store SA
  // params. The contribution from atoms without neighbors depends only on 
  // their vdW radii and P1 parameter so calculate that now.
  SurfInfo SI;
  SI.vdwradii = 0.0;
  SI.P1 = 0.0;
  SI.P2 = 0.0;
  SI.P3 = 0.0;
  SI.P4 = 0.0;
  noNeighborTerm_ = 0.0;
  HeavyAtoms_.clear();
  VDW_.clear();
  SA_Atoms_.clear();
  Params_.clear();
  for (int idx = 0; idx != SoluteMask_.Nselected(); ++idx)
  {
    int atom = SoluteMask_[idx];
    // TODO streamline SurfInfo
    SetAtomLCPO( setup.Top(), atom, &SI );
    mprintf("Atom %i vdw %f\n", atom, SI.vdwradii);
    if (SI.vdwradii > neighborCut_) {
      // Atom has neighbors.
      HeavyAtoms_.push_back( atom );
      VDW_.push_back( SI.vdwradii );
      if ( cmask1.AtomInCharMask( atom ) ) {
        // This is an atom for which we want the SA
        SA_Atoms_.push_back( atom );
        Params_.push_back( SI );
      }
    } else {
      // Atom has no neighbors
      if ( cmask1.AtomInCharMask( atom ) ) {
        // Calculate surface area of atom i
        double vdwi2 = SI.vdwradii * SI.vdwradii;
        double Si = vdwi2 * Constants::FOURPI; 
        mprintf("DBG: AtomNoNbr %i P1 %g Si %g\n", atom, SI.P1, Si);
        noNeighborTerm_ += (SI.P1 * Si);
      }
    }
  } // END loop over solute atoms
  mprintf("\t%zu atoms with neighbors.\n", HeavyAtoms_.size());
  mprintf("\tCalculating SA for %zu atoms with neighbors.\n", SA_Atoms_.size());
  mprintf("\tContribution from atoms with no neighbors is %g\n", noNeighborTerm_);
# ifdef _OPENMP
  // Each thread needs temp. space to store neighbor list and distances for
  // each atom to avoid memory clashes.
# pragma omp parallel
  {
# pragma omp master
  {
  Ineighbor_.resize( omp_get_num_threads() );
  DIJ_.resize( omp_get_num_threads() );
  }
  }
# endif

  return Action::OK;  
}

// Action_Surf::DoAction()
/** Calculate surface area. */
Action::RetType Action_Surf::DoAction(int frameNum, ActionFrame& frm) {
  int idx; // Index into SA_Atoms_ and Params_
  double SA = 0.0;
  int maxIdx = (int)SA_Atoms_.size();
# ifdef _OPENMP
# pragma omp parallel private(idx) reduction(+:SA)
  {
  int mythread = omp_get_thread_num();
  Iarray& ineighbor = Ineighbor_[mythread];
  Darray& Distances_i_j = DIJ_[mythread];
# pragma omp for
# else
  Iarray& ineighbor = Ineighbor_;
  Darray& Distances_i_j = DIJ_;
# endif
  for (idx = 0; idx < maxIdx; idx++)
  {
    int atomi = SA_Atoms_[idx];
    double vdwi = Params_[idx].vdwradii;
    // Search all heavy solute atoms for neighbors of atomi.
    ineighbor.clear();
    Distances_i_j.clear();
    for (unsigned int jdx = 0; jdx != HeavyAtoms_.size(); jdx++)
    {
      int atomj = HeavyAtoms_[jdx];
      if (atomi != atomj) {
        double dij = sqrt( DIST2_NoImage(frm.Frm().XYZ(atomi), frm.Frm().XYZ(atomj)) );
        // Count atoms as neighbors if their vdW radii touch
        if ( (vdwi + VDW_[jdx]) > dij ) {
          ineighbor.push_back( jdx );
          Distances_i_j.push_back( dij );
        }
//        mprintf("SURF_NEIG:  %i %i %f %f %f\n",atomi,atomj,dij,vdwi,VDW_[jdx]);
      }
    }
    // DEBUG - print neighbor list
//    mprintf("SURF: Neighbors for atom %i:",atomi);
//    for (Iarray::const_iterator jt = ineighbor.begin(); jt != ineighbor.end(); jt++)
//      mprintf(" %i",HeavyAtoms_[*jt]);
//    mprintf("\n");
    // Calculate surface area Ai for atomi:
    // Ai = P1*S1 + P2*Sum(Aij) + P3*Sum(Ajk) + P4*Sum(Aij * Sum(Ajk))
    double sumaij = 0.0;
    double sumajk = 0.0;
    double sumaijajk = 0.0;
    double vdwi2 = vdwi * vdwi;
    double Si = vdwi2 * Constants::FOURPI;
    // Loop over all neighbors of atomi (j)
    for (unsigned int mm = 0; mm < ineighbor.size(); mm++)
    {
      double dij = Distances_i_j[mm];
      int jdx = ineighbor[mm];
      int atomj = HeavyAtoms_[jdx];
//      mprintf("i,j %i %i\n",atomi + 1,atomj+1);
      double vdwj = VDW_[jdx];
      double vdwj2 = vdwj * vdwj;
      double tmpaij = vdwi - (dij * 0.5) - ( (vdwi2 - vdwj2)/(2.0 * dij) );
      double aij = Constants::TWOPI * vdwi * tmpaij;
      sumaij += aij;
      // Find which neighbors of atom i (j and k) are themselves neighbors.
      double sumajk_2 = 0.0;
      // FIXME should this start at mm+1?
      for (unsigned int nn = 0; nn < ineighbor.size(); nn++)
      {
        int kdx = ineighbor[nn];
        if (mm != nn)
        {
          int atomk = HeavyAtoms_[kdx];
//          mprintf("i,j,k %i %i %i\n",atomi + 1,atomj+1,atomk+1);
          double vdwk = VDW_[kdx];
          double djk = sqrt(DIST2_NoImage(frm.Frm().XYZ(atomj), frm.Frm().XYZ(atomk)));
//          mprintf("%4s%6i%6i%12.8lf\n","DJK ",atomj+1,atomk+1,djk);
//          mprintf("%6s%6.2lf%6.2lf\n","AVD ",vdwj,vdwk);
          if ( (vdwj + vdwk) > djk ) {
            double vdw2dif = vdwj2 - (vdwk * vdwk);
            double tmpajk = (2.0*vdwj) - djk - (vdw2dif / djk);
            double ajk = Constants::PI*vdwj*tmpajk;
            //tmpajk = vdwj - (djk *0.5) - ( (vdwj2 - (vdwk * vdwk))/(2.0 * djk) );
            //ajk = 2.0 * PI * vdwi * tmpajk;
            //printf("%4s%6i%6i%12.8lf%12.8lf%12.8lf\n","AJK ",(*jt)+1,(*kt)+1,ajk,vdw2dif,tmpajk);
            sumajk += ajk;
            sumajk_2 += ajk;
          }
        }
      } // END loop over neighbor-neighbor pairs of atom i (j to k)
      sumaijajk += (aij * sumajk_2);
      // DEBUG
//      mprintf("%4s%20.8lf %20.8lf %20.8lf\n","AJK ",aij,sumajk,sumaijajk);
    } // END Loop over neighbors of atom i (j)
    SA += ( (Params_[idx].P1 * Si       ) +
            (Params_[idx].P2 * sumaij   ) +
            (Params_[idx].P3 * sumajk   ) +
            (Params_[idx].P4 * sumaijajk) +
            noNeighborTerm_
          ) ;
  } // END loop over atoms in mask (i)
# ifdef _OPENMP
  } // END pragma omp parallel
# endif

  surf_->Add(frameNum, &SA);
    
  return Action::OK;
} 

// -----------------------------------------------------------------------------
// Action_Surf::AssignLCPO()
/** Assign parameters for LCPO method. All radii are incremented by 1.4 Ang. */
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
static void WarnLCPO(NameType const& atype, int atom, int numBonds) {
  mprintf("Warning: Unusual number of bonds for atom %i (%i), type %s.\n",
          atom, numBonds, *atype);
  mprintf("Using default atom parameters.\n");
}

// Action_Surf::SetAtomLCPO()
/** Set up parameters only used in surface area calcs.
  * Adapted from gbsa=1 method in SANDER, mdread.F90
  * \param currentParm The Topology containing atom information. 
  * \param atidx The atom number to set up parameters for.
  * \param SIptr Address to store the SI parameters.
  */
void Action_Surf::SetAtomLCPO(Topology const& currentParm, int atidx, SurfInfo* SIptr) 
{
  const Atom& atom = currentParm[atidx];
  // Get the number of non-H bonded neighbors to this atom
  int numBonds = 0;
  for (Atom::bond_iterator batom = atom.bondbegin(); batom != atom.bondend(); batom++)
    if ( currentParm[ *batom ].Element() != Atom::HYDROGEN )
      ++numBonds;
  char atype0 = toupper(atom.Type()[0]);
  char atype1 = toupper(atom.Type()[1]);
  // TODO: Only set parameters for solute atoms?
  // Set vdw radii and LCPO parameters for this atom
  switch (atom.Element()) {
    case Atom::CARBON:
      if (atom.Nbonds() == 4) {
        switch ( numBonds ) {
          case 1: AssignLCPO(SIptr, 1.70, 0.77887, -0.28063, -0.0012968, 0.00039328); break;
          case 2: AssignLCPO(SIptr, 1.70, 0.56482, -0.19608, -0.0010219, 0.0002658);  break;
          case 3: AssignLCPO(SIptr, 1.70, 0.23348, -0.072627, -0.00020079, 0.00007967); break;
          case 4: AssignLCPO(SIptr, 1.70, 0.00000, 0.00000, 0.00000, 0.00000); break;
          default: WarnLCPO(atom.Type(),atidx + 1,numBonds);
                   AssignLCPO(SIptr, 1.70, 0.77887, -0.28063, -0.0012968, 0.00039328);
        }
      } else {
        switch ( numBonds ) {
          case 2: AssignLCPO(SIptr, 1.70, 0.51245, -0.15966, -0.00019781, 0.00016392); break;
          case 3: AssignLCPO(SIptr, 1.70, 0.070344, -0.019015, -0.000022009, 0.000016875); break;
          default: WarnLCPO(atom.Type(),atidx + 1,numBonds);
                   AssignLCPO(SIptr, 1.70, 0.77887, -0.28063, -0.0012968, 0.00039328);
        }
      }
      break;
    case Atom::OXYGEN:
      if (atype0=='O' && atype1==' ') 
        AssignLCPO(SIptr, 1.60, 0.68563, -0.1868, -0.00135573, 0.00023743);
      else if (atype0=='O' && atype1=='2')
        AssignLCPO(SIptr, 1.60, 0.88857, -0.33421, -0.0018683, 0.00049372);
      else {
        switch ( numBonds ) {
          case 1: AssignLCPO(SIptr, 1.60, 0.77914, -0.25262, -0.0016056, 0.00035071); break;
          case 2: AssignLCPO(SIptr, 1.60, 0.49392, -0.16038, -0.00015512, 0.00016453); break;
          default: WarnLCPO(atom.Type(),atidx + 1,numBonds);
                   AssignLCPO(SIptr, 1.60, 0.77914, -0.25262, -0.0016056, 0.00035071);
        }
      }
      break;
    case Atom::NITROGEN:
      if (atype0=='N' && atype1=='3') {
        switch ( numBonds ) {
          case 1: AssignLCPO(SIptr, 1.65, 0.078602, -0.29198, -0.0006537, 0.00036247); break;
          case 2: AssignLCPO(SIptr, 1.65, 0.22599, -0.036648, -0.0012297, 0.000080038); break;
          case 3: AssignLCPO(SIptr, 1.65, 0.051481, -0.012603, -0.00032006, 0.000024774); break;
          default: WarnLCPO(atom.Type(),atidx + 1,numBonds);
                   AssignLCPO(SIptr, 1.65, 0.078602, -0.29198, -0.0006537, 0.00036247);
        }
      } else {
        switch ( numBonds ) {
          case 1: AssignLCPO(SIptr, 1.65, 0.73511, -0.22116, -0.00089148, 0.0002523); break;
          case 2: AssignLCPO(SIptr, 1.65, 0.41102, -0.12254, -0.000075448, 0.00011804); break;
          case 3: AssignLCPO(SIptr, 1.65, 0.062577, -0.017874, -0.00008312, 0.000019849); break;
          default: WarnLCPO(atom.Type(),atidx + 1,numBonds);
                   AssignLCPO(SIptr, 1.65, 0.078602, -0.29198, -0.0006537, 0.00036247);
        }
      }
      break;
    case Atom::SULFUR:
      if (atype0=='S' && atype1=='H') 
        AssignLCPO(SIptr, 1.90, 0.7722, -0.26393, 0.0010629, 0.0002179);
      else 
        AssignLCPO(SIptr, 1.90, 0.54581, -0.19477, -0.0012873, 0.00029247);
      break;
    case Atom::PHOSPHORUS:
      switch ( numBonds ) {
        case 3: AssignLCPO(SIptr, 1.90, 0.3865, -0.18249, -0.0036598, 0.0004264); break;
        case 4: AssignLCPO(SIptr, 1.90, 0.03873, -0.0089339, 0.0000083582, 0.0000030381); break;
        default: WarnLCPO(atom.Type(),atidx + 1,numBonds);
                 AssignLCPO(SIptr, 1.90, 0.3865, -0.18249, -0.0036598, 0.0004264);
      }
      break;
    case Atom::HYDROGEN:
      AssignLCPO(SIptr, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000); 
      break;
    default:
      if (atype0=='Z') 
        AssignLCPO(SIptr, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000);
      else if (atype0=='M' && atype1=='G') 
        //  Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
        //  Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
        //  Mg radius = 1.45A: Aqvist 1992
        //  The P1-4 values were taken from O.sp3 with two bonded 
        //  neighbors -> O has the smallest van der Waals radius 
        //  compared to all other elements which had been parametrized
        AssignLCPO(SIptr, 1.18, 0.49392, -0.16038, -0.00015512, 0.00016453);
      else if (atype0=='F')
        AssignLCPO(SIptr, 1.47, 0.68563, -0.1868, -0.00135573, 0.00023743);
      else {
        mprintf("Warning: Using carbon SA parms for unknown atom %i type %s\n",
                atidx + 1, *(atom.Type()));
        AssignLCPO(SIptr, 1.70, 0.51245, -0.15966, -0.00019781, 0.00016392);
      }
  } // END switch atom.Element() 
}

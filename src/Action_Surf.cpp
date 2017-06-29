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
Action_Surf::Action_Surf() : surf_(0) {} 

// Action_Surf::Help()
void Action_Surf::Help() const {
  mprintf("\t<name> [<mask1>] [out <filename>]\n"
          "  Calculate LCPO surface area of atoms in <mask1>, or all solute\n"
          "  atoms if no mask specified.\n");
}

// Action_Surf::Init()
Action::RetType Action_Surf::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get keywords
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs);

  // Get Masks
  std::string maskexp = actionArgs.GetMaskNext();
  if (!maskexp.empty()) Mask1_.SetMaskString( maskexp );

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
  SoluteAtoms_.ResetMask();
  SoluteAtoms_.SetNatoms( setup.Top().Natom() );
  if ( setup.Top().Nmol() > 0) {
    for (int at = 0; at != setup.Top().Natom(); at++)
    {
      Molecule const& mol = setup.Top().Mol( setup.Top()[at].MolNum() );
      if (!mol.IsSolvent() && mol.NumAtoms() > 1)
        SoluteAtoms_.AddSelectedAtom( at );
    }
  } else {
    mprintf("Warning: No molecule info in '%s'. Adding all atoms.\n");
    for (int at = 0; at != setup.Top().Natom(); at++)
      SoluteAtoms_.AddSelectedAtom( at );
  }
  if (!Mask1_.MaskStringSet())
    Mask1_ = SoluteAtoms_;
  mprintf("\t%i solute atoms. Calculating LCPO surface area for %i atoms.\n",
          SoluteAtoms_.Nselected(), Mask1_.Nselected());

  return Action::OK;  
}

// Action_Surf::DoAction()
/** Calculate surface area. */
Action::RetType Action_Surf::DoAction(int frameNum, ActionFrame& frm) {

  return Action::ERR;
} 
/*
// -----------------------------------------------------------------------------
// Action_Surf::AssignLCPO()
* Assign parameters for LCPO method. All radii are incremented by 1.4 Ang.

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
* Set up parameters only used in surface area calcs.
  * Adapted from gbsa=1 method in SANDER, mdread.F90
  * \param currentParm The Topology containing atom information. 
  * \param atidx The atom number to set up parameters for.
  * \param SIptr Address to store the SI parameters.

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
*/

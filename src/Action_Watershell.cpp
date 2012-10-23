#include "Action_Watershell.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Watershell::Action_Watershell() :
  lowerCutoff_(0),
  upperCutoff_(0),
  CurrentParm_(0)
{ }

void Action_Watershell::Help() {
  mprintf("watershell <solutemask> <filename> [lower <lower cut>] [upper <upper cut>]\n");
  mprintf("           [noimage] [<solventmask>]\n");
  mprintf("\tCalculate # of waters in 1st and 2nd solvation shells.\n");
}

// Action_Watershell::init()
Action::RetType Action_Watershell::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  InitImaging( !actionArgs.hasKey("noimage") );

  ArgList::ConstArg maskexpr = actionArgs.getNextMask();
  if (maskexpr == NULL) {
    mprinterr("Error: WATERSHELL: Solute mask must be specified.\n");
    return Action::ERR;
  }
  soluteMask_.SetMaskString( maskexpr );

  filename_ = actionArgs.GetStringNext();
  if (filename_.empty()) {
    mprinterr("Error: WATERSHELL: Output filename must be specified.\n");
    return Action::ERR;
  }

  lowerCutoff_ = actionArgs.getKeyDouble("lower", 3.4);
  upperCutoff_ = actionArgs.getKeyDouble("upper", 5.0);

  // Check for solvent mask
  solventmaskexpr_ = actionArgs.GetMaskNext();

  mprintf("    WATERSHELL: Output to %s\n",filename_.c_str());
  if (!UseImage())
    mprintf("                Imaging is disabled.\n");
  mprintf("                The first shell will contain water < %.3lf angstroms from\n",
          lowerCutoff_);
  mprintf("                the solute; the second shell < %.3lf angstroms...\n",
          upperCutoff_);
  mprintf("                Solute atoms will be specified by [%s]\n",soluteMask_.MaskString());
  if (!solventmaskexpr_.empty()) {
    mprintf("                Solvent atoms will be specified by [%s]\n",
            solventmaskexpr_.c_str());
    solventMask_.SetMaskString( solventmaskexpr_ );
  }

  // Pre-square upper and lower cutoffs
  lowerCutoff_ *= lowerCutoff_;
  upperCutoff_ *= upperCutoff_;

  return Action::OK;
}

// Action_Watershell::setup()
/** Set up solute and solvent masks. If no solvent mask was specified use 
  * solvent information in the current topology.
  */
Action::RetType Action_Watershell::Setup(Topology* currentParm, Topology** parmAddress) {
  // Set up solute mask
  if (currentParm->SetupIntegerMask( soluteMask_ )) return Action::ERR;
  if ( soluteMask_.None() ) {
    mprinterr("Error: WATERSHELL: No atoms in solute mask [%s].\n",soluteMask_.MaskString());
    return Action::ERR;
  }
  // Set up solvent mask
  if (!solventmaskexpr_.empty()) {
    if (currentParm->SetupIntegerMask( solventMask_ )) return Action::ERR;
  } else {
    solventMask_.ResetMask();
    for (Topology::mol_iterator mol = currentParm->MolStart();
                                mol != currentParm->MolEnd(); ++mol)
    {
      if ( (*mol).IsSolvent() )
        solventMask_.AddAtomRange( (*mol).BeginAtom(), (*mol).EndAtom() );
    }
  }
  if ( solventMask_.None() ) {
    if (!solventmaskexpr_.empty())
      mprinterr("Error: WATERSHELL: No solvent atoms selected by mask [%s]\n",
                solventmaskexpr_.c_str());
    else
      mprinterr("Error: WATERSHELL: No solvent atoms in topology %s\n",currentParm->c_str());
    return Action::ERR;
  }
  SetupImaging( currentParm->BoxType() );
  // Create space for residues
  activeResidues_.resize( currentParm->Nres(), 0 );
  // Store current Parm
  CurrentParm_ = currentParm;
  return Action::OK;    
}

// Action_Watershell::action()
Action::RetType Action_Watershell::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  double ucell[9], recip[9];
  Vec3 boxXYZ(currentFrame->BoxX(), currentFrame->BoxY(), currentFrame->BoxZ() );
 
  if (ImageType()==NONORTHO) currentFrame->BoxToRecip(ucell,recip);

  // Loop over solute atoms
  for (AtomMask::const_iterator solute_at = soluteMask_.begin();
                                solute_at != soluteMask_.end(); ++solute_at)
  {
    // Loop over solvent atoms
    for (AtomMask::const_iterator solvent_at = solventMask_.begin();
                                  solvent_at != solventMask_.end(); ++solvent_at)
    {
      // Figure out which solvent residue this is
      int currentRes = (*CurrentParm_)[ *solvent_at].ResNum();
      // If residue is not yet marked as 1st shell, calc distance
      if ( activeResidues_[currentRes] < 2 ) {
        double dist = DIST2(currentFrame->XYZ(*solute_at), currentFrame->XYZ(*solvent_at), 
                            ImageType(), boxXYZ, ucell, recip );
        // Less than upper, 2nd shell
        if (dist < upperCutoff_) {
          activeResidues_[currentRes] = 1;
          // Less than lower, 1st shell
          if (dist < lowerCutoff_) 
            activeResidues_[currentRes] = 2;
        }
      }
    } // END loop over solvent atoms
  } // END loop over solute atoms

  // Now each residue is marked 0 (no shell), 1 (second shell), 2 (first shell)
  // TODO: Use datasets instead
  lower_.push_back( 0 );
  upper_.push_back( 0 );
  for (std::vector<int>::iterator shell = activeResidues_.begin();
                                  shell != activeResidues_.end(); ++shell)
  {
    if ( *shell > 0 ) {
      ++(upper_.back());
      if ( *shell > 1 ) {
        ++(lower_.back());
      }
    }
    // Reset for next pass
    *shell = 0;
  }

  return Action::OK;
}

// Action_Watershell::print()
void Action_Watershell::Print() {
  CpptrajFile outfile;
  if (outfile.OpenWrite( filename_ )) return;

  unsigned int frame = 1;
  std::vector<int>::iterator L = lower_.begin();
  for (std::vector<int>::iterator U = upper_.begin(); U != upper_.end(); ++U) { 
    outfile.Printf("%i %i %i\n", frame++, *L, *U);
    ++L;
  }
  
  outfile.CloseFile();
}


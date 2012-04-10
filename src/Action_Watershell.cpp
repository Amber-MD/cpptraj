#include "Action_Watershell.h"
#include "CpptrajStdio.h"

Watershell::Watershell() :
  solventmaskexpr_(0),
  //visits_(0),
  lowerCutoff_(0),
  upperCutoff_(0),
  filename_(0)
{ }

/** Expected call: 
  * watershell <solutemask> <filename> [lower <lower cut>] [upper <upper cut>] 
  *            [noimage] [<solventmask>]
  *      
  */
int Watershell::init() {
  useImage = !actionArgs.hasKey("noimage");

  char *maskexpr = actionArgs.getNextMask();
  if (maskexpr == NULL) {
    mprinterr("Error: WATERSHELL: Solute mask must be specified.\n");
    return 1;
  }
  soluteMask_.SetMaskString( maskexpr );

  filename_ = actionArgs.getNextString();
  if (filename_ == NULL) {
    mprinterr("Error: WATERSHELL: Output filename must be specified.\n");
    return 1;
  }

  lowerCutoff_ = actionArgs.getKeyDouble("lower", 3.4);
  upperCutoff_ = actionArgs.getKeyDouble("upper", 5.0);

  // Check for solvent mask
  solventmaskexpr_ = actionArgs.getNextMask();

  mprintf("    WATERSHELL: Output to %s\n",filename_);
  if (!useImage)
    mprintf("                Imaging is disabled.\n");
  mprintf("                The first shell will contain water < %.3lf angstroms from\n",
          lowerCutoff_);
  mprintf("                the solute; the second shell < %.3lf angstroms...\n",
          upperCutoff_);
  mprintf("                Solute atoms will be specified by [%s]\n",soluteMask_.MaskString());
  if (solventmaskexpr_!=NULL) {
    mprintf("                Solvent atoms will be specified by [%s]\n",solventmaskexpr_);
    solventMask_.SetMaskString( solventmaskexpr_ );
  }

  // Pre-square upper and lower cutoffs
  lowerCutoff_ *= lowerCutoff_;
  upperCutoff_ *= upperCutoff_;

  return 0;
}

int Watershell::setup() {
  // Set up solute mask
  if (currentParm->SetupIntegerMask( soluteMask_ )) return 1;
  if ( soluteMask_.None() ) {
    mprinterr("Error: WATERSHELL: No atoms in solute mask [%s].\n",soluteMask_.MaskString());
    return 1;
  }
  // Set up solvent mask
  if (solventmaskexpr_ != NULL) {
    if (currentParm->SetupIntegerMask( solventMask_ )) return 1;
  } else {
    solventMask_.ResetMask();
    for (Topology::mol_iterator mol = currentParm->SolventStart();
                                mol != currentParm->SolventEnd(); ++mol)
    {
      solventMask_.AddAtomRange( (*mol).BeginAtom(), (*mol).EndAtom() );
    }
  }
  if ( solventMask_.None() ) {
    if (solventmaskexpr_!=NULL)
      mprinterr("Error: WATERSHELL: No solvent atoms selected by mask [%s]\n",solventmaskexpr_);
    else
      mprinterr("Error: WATERSHELL: No solvent atoms in topology %s\n",currentParm->c_str());
    return 1;
  }
  // Create space for residues
  activeResidues_.resize( currentParm->Nres(), 0 );
  return 0;    
}

int Watershell::action() {
  double ucell[9], recip[9];

  if (imageType==NONORTHO) currentFrame->BoxToRecip(ucell,recip);

  // Loop over solute atoms
  for (AtomMask::const_iterator solute_at = soluteMask_.begin();
                                solute_at != soluteMask_.end(); ++solute_at)
  {
    // Loop over solvent atoms
    for (AtomMask::const_iterator solvent_at = solventMask_.begin();
                                  solvent_at != solventMask_.end(); ++solvent_at)
    {
      // Figure out which solvent residue this is
      int currentRes = (*currentParm)[ *solvent_at].ResNum();
      // If residue is not yet marked as 1st shell, calc distance
      if ( activeResidues_[currentRes] < 2 ) {
        double dist = currentFrame->DIST2( *solute_at, *solvent_at, imageType,
                                           ucell, recip );
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
  }

  return 0;
}

void Watershell::print() {
  CpptrajFile outfile;
  if (outfile.SetupWrite( filename_, debug )) return;
  outfile.OpenFile();

  unsigned int frame = 1;
  std::vector<int>::iterator L = lower_.begin();
  for (std::vector<int>::iterator U = upper_.begin(); U != upper_.end(); ++U) 
    outfile.Printf("%i %i %i\n", frame++, *L, *U);
  
  outfile.CloseFile();
}


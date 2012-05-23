#include "Action_Unwrap.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Unwrap::Action_Unwrap() :
  RefParm_(0),
  orthogonal_(false)
{ }

int Action_Unwrap::init() {
  // Get reference
  int refindex = -1;
  if (actionArgs.hasKey("reference"))
    refindex = 0;
  else {
    refindex = actionArgs.getKeyInt("refindex", -1);
    std::string refname = actionArgs.GetStringKey("ref");
    if (!refname.empty()) {
      refindex = FL->FindName( refname );
      if (refindex == -1) {
        mprinterr("Error: unwrap: Reference [%s] not found.\n", refname.c_str());
        return 1;
      }
    }
  }
  // If refindex is not -1, attempt to get reference frame.
  if (refindex > -1) {
    Frame *tempframe = FL->GetFrame(refindex);
    if (tempframe==NULL) {
      mprinterr("Error: unwrap: Could not get reference index %i\n", refindex);
      return 1;
    }
    RefFrame_ = *tempframe;
    // Get reference parm for frame
    RefParm_ = FL->GetFrameParm(refindex);
    if (RefParm_ == NULL) {
      mprinterr("Error: unwrap: Could not get reference parm for frame %s\n",
                FL->FrameName(refindex));
      return 1;
    }
  }

  // Get mask string
  char* maskexpr = actionArgs.getNextMask();
  mask_.SetMaskString( maskexpr );

  mprintf("    UNWRAP: (%s), reference is ", mask_.MaskString());
  if ( refindex > -1)
    mprintf("%s", FL->FrameName(refindex));
  else
    mprintf("first frame.");
  mprintf("\n");

  return 0;
}

int Action_Unwrap::setup() {
  // Ensure same number of atoms in current parm and ref parm
  if ( RefParm_!=0 ) {
    if ( currentParm->Natom() != RefParm_->Natom() ) {
      mprinterr("Error: unwrap: # atoms in reference parm %s is not\n", RefParm_->c_str());
      mprinterr("Error:         equal to # atoms in parm %s\n", currentParm->c_str());
      return 1;
    }
  }

  // Setup mask
  if ( currentParm->SetupIntegerMask( mask_ ) ) return 1;
  if (mask_.None()) {
    mprinterr("Error: unwrap: No atoms selected.\n");
    return 1;
  }
  mprintf("\t[%s] %i atoms selected.\n", mask_.MaskString(), mask_.Nselected());

  // Check box type
  if (currentParm->BoxType()==Box::NOBOX) {
    mprintf("Error: unwrap: Parm %s does not contain box information.\n",
            currentParm->c_str());
    return 1;
  }
  orthogonal_ = false;
  if (currentParm->BoxType()==Box::ORTHO)
    orthogonal_ = true;

  return 0;
}

int Action_Unwrap::action() {
  // Set reference structure if not already set
  if (RefParm_ == 0) {
    RefParm_ = currentParm;
    RefFrame_ = *currentFrame;
    return 0;
  }
 
  if (orthogonal_)
    currentFrame->UnwrapOrtho( RefFrame_, mask_ );
  else 
    currentFrame->UnwrapNonortho( RefFrame_, mask_ );

  return 0;
} 

#include "Action_Unwrap.h"
#include "CpptrajStdio.h"
#include "ImageRoutines.h"

// CONSTRUCTOR
Action_Unwrap::Action_Unwrap() :
  RefParm_(0),
  orthogonal_(false)
{ }

void Action_Unwrap::Help() {

}

Action::RetType Action_Unwrap::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
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
        return Action::ERR;
      }
    }
  }
  // If refindex is not -1, attempt to get reference frame.
  if (refindex > -1) {
    Frame *tempframe = FL->GetFrame(refindex);
    if (tempframe==NULL) {
      mprinterr("Error: unwrap: Could not get reference index %i\n", refindex);
      return Action::ERR;
    }
    RefFrame_ = *tempframe;
    // Get reference parm for frame
    RefParm_ = FL->GetFrameParm(refindex);
    if (RefParm_ == NULL) {
      mprinterr("Error: unwrap: Could not get reference parm for frame %s\n",
                FL->FrameName(refindex));
      return Action::ERR;
    }
  }

  // Get mask string
  mask_.SetMaskString( actionArgs.getNextMask() );

  mprintf("    UNWRAP: (%s), reference is ", mask_.MaskString());
  if ( refindex > -1)
    mprintf("%s", FL->FrameName(refindex));
  else
    mprintf("first frame.");
  mprintf("\n");

  return Action::OK;
}

Action::RetType Action_Unwrap::Setup(Topology* currentParm, Topology** parmAddress) {
  // Ensure same number of atoms in current parm and ref parm
  if ( RefParm_!=0 ) {
    if ( currentParm->Natom() != RefParm_->Natom() ) {
      mprinterr("Error: unwrap: # atoms in reference parm %s is not\n", RefParm_->c_str());
      mprinterr("Error:         equal to # atoms in parm %s\n", currentParm->c_str());
      return Action::ERR;
    }
  }

  // Setup mask
  if ( currentParm->SetupIntegerMask( mask_ ) ) return Action::ERR;
  if (mask_.None()) {
    mprinterr("Error: unwrap: No atoms selected.\n");
    return Action::ERR;
  }
  mprintf("\t[%s] %i atoms selected.\n", mask_.MaskString(), mask_.Nselected());

  // Check box type
  if (currentParm->BoxType()==Box::NOBOX) {
    mprintf("Error: unwrap: Parm %s does not contain box information.\n",
            currentParm->c_str());
    return Action::ERR;
  }
  orthogonal_ = false;
  if (currentParm->BoxType()==Box::ORTHO)
    orthogonal_ = true;
  // Use current parm as reference if not already set
  if (RefParm_ == 0)
    RefParm_ = currentParm;
  return Action::OK;
}

Action::RetType Action_Unwrap::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  double ucell[9], recip[9];
  // Set reference structure if not already set
  if (RefFrame_.empty()) {
    RefFrame_ = *currentFrame;
    return Action::OK;
  }
 
  if (orthogonal_)
    UnwrapOrtho( *currentFrame, RefFrame_, mask_ );
  else {
    currentFrame->BoxToRecip( ucell, recip );
    UnwrapNonortho( *currentFrame, RefFrame_, mask_, ucell, recip );
  }

  return Action::OK;
} 

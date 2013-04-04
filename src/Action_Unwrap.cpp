#include "Action_Unwrap.h"
#include "CpptrajStdio.h"
#include "ImageRoutines.h"

// CONSTRUCTOR
Action_Unwrap::Action_Unwrap() :
  RefParm_(0),
  orthogonal_(false)
{ }

void Action_Unwrap::Help() {
  mprintf("\t[{reference | ref <refname> | refindex <#>}] [<mask>]\n");
  mprintf("\tReverse of 'image'; unwrap coordinates in <mask> according\n");
  mprintf("\tto a reference structure.\n");
}

Action::RetType Action_Unwrap::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get reference
  ReferenceFrame REF = FL->GetFrameFromArgs( actionArgs );
  if (REF.error()) return Action::ERR;
  if (!REF.empty()) {
    RefFrame_ = *(REF.Coord());
    // Get reference parm for frame
    RefParm_ = REF.Parm();
  }

  // Get mask string
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("    UNWRAP: (%s), reference is ", mask_.MaskString());
  if ( !REF.empty())
    mprintf("%s", REF.FrameName().c_str());
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
  Matrix_3x3 ucell, recip;
  // Set reference structure if not already set
  if (RefFrame_.empty()) {
    RefFrame_ = *currentFrame;
    return Action::OK;
  }
 
  if (orthogonal_)
    UnwrapOrtho( *currentFrame, RefFrame_, mask_ );
  else {
    currentFrame->BoxCrd().ToRecip( ucell, recip );
    UnwrapNonortho( *currentFrame, RefFrame_, mask_, ucell, recip );
  }

  return Action::OK;
} 

#include "ActionReference.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
ActionReference::ActionReference() {
  RefParm_ = NULL;
  fitReference_ = false;
  useRefMass_ = false;
  refMode_ = NOREF;
  allowed_.assign(5, true);
  refParmSetup_ = false;
}

// ActionReference::SetFirst()
/** Set Reference Mode to FIRST, which requires no additional setup.
  */
// NOTE: Necessary?
void ActionReference::SetFirst(bool nofit, char *mask0, bool useMassIn) {
  refMode_ = FIRST;
  fitReference_ = (!nofit);
  useRefMass_ = useMassIn;
  RefMask_.SetMaskString(mask0);
}

// ActionReference::RefInfo()
/** Print information on how reference has been set up. */
void ActionReference::RefInfo() {
  //mprintf("\tReference is");
  if (refMode_==FIRST)
    mprintf(" first frame");
  else if (refMode_==REFTRAJ) 
    mprintf(" trajectory %s",RefTraj_.TrajName());
  else {
    mprintf(" reference frame");
    if (!refname_.empty())
      mprintf(" %s",refname_.c_str());
  }
  mprintf(" (%s)",RefMask_.MaskString());
}

// ActionReference::RefNselected()
int ActionReference::RefNselected() {
  return RefMask_.Nselected();
}

// ActionReference::GetRefParm()
Topology *ActionReference::GetRefParm() {
  return RefParm_;
}

// ActionReference::RefInit()
/** Initialize reference information from argument list. This should be called
  * after all other keywords for the action have been processed. 
  */
int ActionReference::RefInit(bool nofit, bool useMassIn, char *defaultMask,
                             ArgList &actionArgs, FrameList *RefFrameList,
                             TopologyList *parmList, double *RefTrans) 
{
  int refindex = -1;
  char *referenceName = NULL;
  char *reftraj = NULL;
  
  fitReference_ = (!nofit);
  useRefMass_ = useMassIn;
  // If fitReference RefTrans cannot be NULL
  if (fitReference_ && RefTrans==NULL) {
    mprinterr("Internal Error: RefInit: fitReference specified but RefTrans is NULL.\n");
    return 1;
  }
  // Check for all reference keywords, priority FIRST > REF > REFTRAJ
  if (allowed_[FIRST] && actionArgs.hasKey("first"))
    refMode_ = FIRST;
  else {
    // Check for ref keywords if allowed
    if (allowed_[REF]) { 
      referenceName = actionArgs.getKeyString("ref",NULL);
      refindex = actionArgs.getKeyInt("refindex",-1);
      // For compatibility with ptraj, keyword 'reference' == 'refindex 0'
      if (actionArgs.hasKey("reference")) refindex = 0;
    }
    // Check for reftraj keywords if allowed
    if (allowed_[REFTRAJ]) {
      reftraj = actionArgs.getKeyString("reftraj",NULL);
      // Get reference parm for traj
      RefParm_ = parmList->GetParm(actionArgs);
    }
    if (referenceName!=NULL || refindex!=-1) 
      refMode_ = REF;
    else if (reftraj != NULL)
      refMode_ = REFTRAJ;
  } 

  // Get the reference mask string. If no reference mask specified,
  // use the given default mask.
  char *maskRef = actionArgs.getNextMask();
  if (maskRef==NULL) maskRef = defaultMask;
  RefMask_.SetMaskString(maskRef);

  // Check reference structure mode, default to first
  if (refMode_ == NOREF) {
    mprintf("    Warning: No reference structure given. Defaulting to first.\n");
    refMode_ = FIRST;
  }

  // If not using first frame, set up reference now
  if (refMode_ == REF) {
    // Attempt to get the reference index by name/tag
    if (referenceName != NULL)
      refindex = RefFrameList->FindName(referenceName);
    // Get reference frame by index
    // TODO: Convert FrameList to return frame reference?
    Frame *TempFrame = RefFrameList->GetFrame(refindex);
    if (TempFrame==NULL) {
      mprinterr("Error: Could not get reference index %i\n",refindex);
      return 1;
    }
    // Get reference parm for frame
    RefParm_ = RefFrameList->GetFrameParm(refindex);
    if (RefParm_==NULL) {
      mprinterr("Error: Could not get parm for frame %s\n",RefFrameList->FrameName(refindex));
      return 1;
    }
    // Set up reference frame
    if (SetRefMask()) return 1;
    SetRefStructure( *TempFrame, RefTrans, RefFrameList->FrameName(refindex));
  } else if (refMode_ == REFTRAJ) {
     if (RefParm_==NULL) {
        mprinterr("Error: Could not get parm for reftraj %s.\n",reftraj);
        return 1;
      }
      // Set up reference traj
      if (SetRefMask()) return 1;
      if (SetRefTraj(reftraj)) return 1;
  }
  return 0;
}

// ActionReference::RefSetup()
int ActionReference::RefSetup(Topology *currentParm) {
  if (refMode_==FIRST) {
    RefParm_ = currentParm;
    if (SetRefMask()) return 1;
  }
  return 0;
}

// ActionReference::RefAction()
/** Set first frame or get next reference frame in trajectory. */
void ActionReference::RefAction(Frame *currentFrame, double *RefTrans) {
  if (refMode_ == FIRST) {
    SetRefStructure( *currentFrame, RefTrans, NULL);
    refMode_ = FIRST_DONE;
  } else if (refMode_ == REFTRAJ) {
    RefTraj_.GetNextFrame(RefFrame_);
    SelectedRef_.SetCoordinates(RefFrame_, RefMask_);
    if (fitReference_)
      SelectedRef_.CenterReference(RefTrans, useRefMass_);
  }
}

/// Set reference mask based on reference parm. Allocate selected ref frame.
int ActionReference::SetRefMask() {
  if (RefParm_==NULL) return 1;
  // Only allow parm to be set up once for now
  if (refParmSetup_) return 0;
  // Set up reference mask
  if (RefParm_->SetupIntegerMask( RefMask_ )) return 1;
  if (RefMask_.None()) {
    mprinterr("Error: SetRefParm: Parm %s, no atoms in reference mask [%s].\n",
            RefParm_->c_str(), RefMask_.MaskString());
    return 1;
  }
  // Check if reference parm has masses
  // NOTE: Mass is now always set to 1 if not read in so this only depends
  //       on what the action set useMass to.
  /*if (useRefMass_ && RefParm_->mass==NULL) {
    mprintf("    Warning: SetRefParm: Ref Parmtop %s does not contain mass info.\n",
            RefParm_->c_str());
    mprintf("             Geometric center will be used instead.\n");
    useRefMass_=false;
  }*/
  // Allocate frame for selected reference atoms
  SelectedRef_.SetupFrameFromMask(RefMask_, RefParm_->Mass());
  //mprintf("DEBUG: RefMask has %i atoms\n",RefMask.Nselected);
  refParmSetup_=true;
  return 0;
}

/// Set selected ref coordinates from frameIn based on RefMask.
void ActionReference::SetRefStructure(Frame &frameIn,double *Trans,const char* nameIn) 
{
  RefFrame_ = frameIn;
  SelectedRef_.SetCoordinates(RefFrame_, RefMask_);
  if (fitReference_)
    SelectedRef_.CenterReference(Trans, useRefMass_);
  if (nameIn!=NULL)
    refname_.assign(nameIn);
  else
    refname_.clear();
}

/// Set reference trajectory, open.
int ActionReference::SetRefTraj(char *filenameIn) {
  if (filenameIn==NULL) return 1;
  if (RefParm_==NULL) return 1;
  if (RefTraj_.SetupRead(filenameIn,NULL,RefParm_)) {
    mprinterr("Error: SetRefTraj: Could not set up reftraj %s.\n",filenameIn);
    return 1;
  }
  RefFrame_.SetupFrameV(RefParm_->Natom(), RefParm_->Mass(), RefTraj_.HasVelocity());
  if (RefTraj_.BeginTraj(false)) {
    mprinterr("Error: Rmsd: Could not open reference trajectory.\n");
    return 1;
  }
  return 0;
}


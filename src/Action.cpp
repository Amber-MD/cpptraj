#include "Action.h"

// CONSTRUCTOR
Action::Action() {
  //currentType=UNKNOWN_ACTION; 
  A=NULL; 
  currentFrame=0; 
  noInit=0; 
  noSetup=0; 
  P=NULL;
  F=NULL;
  debug=0;
}

// DESTRUCTOR - virtual since this class is inherited
Action::~Action() {
  //fprintf(stderr,"Action Destructor.\n"); 
  if (A!=NULL) delete A;
}

/*
 * Action::setArg()
 * Set the argument list
 */
void Action::setArg(ArgList *inA) { 
  A=inA->Copy();       
}

/*
 * Action::ResetArg()
 * Reset arguments in the argument list
 */
void Action::ResetArg() { 
  A->Reset();          
}

/*
 * Action::Name()
 *  Print the command that calls the action
 */
char *Action::Name() { 
  return A->Command(); 
}

/*
 * Action::CmdLine()
 * Print the command and all args
 */
char *Action::CmdLine() { 
  return A->ArgLine(); 
}

/*
 * Action::Init()
 * Allocate non-parm dependent things like I/O files and datasets. Sets
 * the pointers to DataSetList, DataFileList and reference FrameList from
 * PtrajState. Called before trajectory processing. Calls the actions 
 * internal init() routine.
 */
int Action::Init(DataSetList *DSLin, FrameList *FLin, DataFileList *DFLin, 
                 ParmFileList *PFLin, int debugIn) {
  int err;

  DSL=DSLin;
  FL=FLin;
  DFL=DFLin;
  PFL=PFLin;
  debug=debugIn;
  err = this->init();
  // Check for unhandled keywords
  this->A->CheckForMoreArgs();

  return ( err );
}

/*
 * Action::Setup()
 * Allocate parm-dependent things like masks. Sets the pointer to the current
 * parmtop and allows it to be modified by the action. Called every time traj 
 * changes. Calls the actions internal setup() routine.
 */
int Action::Setup(AmberParm **ParmAddress) {
  P = *ParmAddress;
  if (this->setup()) return 1;
  // Set the value of parm address in case parm was changed, e.g. in strip
  *ParmAddress = P;
  return 0;
}

/*
 * Action::DoAction() 
 * Perform the action. Called every time a frame is to be processed. Sets the
 * pointer to the current frame and allows it to be modified by the action.
 * Calls the actions internal action() routine.
 */
int Action::DoAction(Frame **FrameAddress, int frameIn) {
  F = *FrameAddress;
  currentFrame = frameIn;
  if (this->action()) return 1;
  // Set the value of frame address in case frame was changed, e.g. in strip
  *FrameAddress = F;
  return 0;
}


#include "Action.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action::Action() {
  P=NULL;
  F=NULL;
  DSL=NULL;
  DFL=NULL;
  PFL=NULL;
  FL=NULL;
  useMass=false;
  debug=0;
  currentFrame=0; 
  noInit=0; 
  noSetup=0; 
}

// DESTRUCTOR
Action::~Action() {
  //fprintf(stderr,"Action Destructor.\n"); 
}

/* Action::SetArg()
 * Set the argument list
 */
void Action::SetArg( const ArgList &inA) { 
  actionArgs = inA;
}

/* Action::ActionCommand()
 *  Print the command that calls the action
 */
const char *Action::ActionCommand() { 
  return actionArgs.Command(); 
}

/* Action::CmdLine()
 * Print the command and all args
 */
const char *Action::CmdLine() { 
  return actionArgs.ArgLine(); 
}

/* Action::Init()
 * Allocate non-parm dependent things like I/O files and datasets. Sets
 * the pointers to DataSetList, DataFileList and reference FrameList from
 * CpptrajState. Called before trajectory processing. Calls the actions 
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
  actionArgs.CheckForMoreArgs();

  return ( err );
}

/* Action::Setup()
 * Allocate parm-dependent things like masks. Sets the pointer to the current
 * parmtop and allows it to be modified by the action. Called every time traj 
 * changes. Calls the actions internal setup() routine.
 */
int Action::Setup(AmberParm **ParmAddress) {
  int err;
  
  P = *ParmAddress;
  // If useMass, check that parm actually has masses.
  if (P->mass==NULL && useMass) {
    mprintf("    Warning: %s: Mass for this parm is NULL.\n",actionArgs.Command());
    mprintf("             Geometric center will be used instead of center of mass.\n");
    useMass=false;
  }
  err = this->setup();
  if (err) return err;
  // Set the value of parm address in case parm was changed, e.g. in strip
  *ParmAddress = P;
  return 0;
}

/* Action::DoAction() 
 * Perform the action. Called every time a frame is to be processed. Sets the
 * pointer to the current frame and allows it to be modified by the action.
 * Calls the actions internal action() routine.
 */
int Action::DoAction(Frame **FrameAddress, int frameIn) {
  int err;

  F = *FrameAddress;
  currentFrame = frameIn;
  err = this->action();
  if (err) return err;
  // Set the value of frame address in case frame was changed, e.g. in strip
  *FrameAddress = F;
  return 0;
}


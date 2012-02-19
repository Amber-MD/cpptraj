#include "Action.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action::Action() {
  isSeparate=false;
  currentParm=NULL;
  currentFrame=NULL;
  DSL=NULL;
  DFL=NULL;
  PFL=NULL;
  FL=NULL;
  activeRefFrame=NULL;
  activeReference=NULL;
  activeReferenceOriginalValue=NULL;
  useMass=false;
  useMassOriginalValue=false;
  useImage=false;
  useImageOriginalValue=false;
  imageType = NOBOX;
  debug=0;
  frameNum=0; 
  noInit=false; 
  noSetup=false; 
}

// DESTRUCTOR
Action::~Action() {
  //fprintf(stderr,"Action Destructor.\n"); 
}

// Action::SetArg()
/** \param inA argument list for the action
  */
void Action::SetArg( const ArgList &inA) { 
  actionArgs = inA;
}

// Action::ActionCommand()
const char *Action::ActionCommand() { 
  return actionArgs.Command(); 
}

// Action::CmdLine()
const char *Action::CmdLine() { 
  return actionArgs.ArgLine(); 
}

// Action::Init()
/** Get all relevant arguments from the
  * argument list, set up any datasets which will be available for 
  * analysis after trajectory processing, and set memory references 
  * to the master dataset, datafile, reference frame, and parm file lists.
  * Actions can have datasets not related to master datasets, e.g.
  * in the NAstruct action, since data is stored for each nucleobase at
  * each frame its output would not match up with other actions in the 
  * master dataset list, so it has its own dataset list.
  * \param DSLin pointer to the master DataSetList
  * \param FLin pointer to the master Reference FrameList
  * \param DFLin pointer to the master DataFileList
  * \param PFLin pointer to the master ParmFileList
  * \param debugIn Debug level that action should be set to
  */
int Action::Init(DataSetList *DSLin, FrameList *FLin, DataFileList *DFLin, 
                 ParmFileList *PFLin, int debugIn) {
  int err;

  DSL=DSLin;
  FL=FLin;
  // Set the active reference frame and coords
  activeRefFrame = FL->ActiveReference();
  if (activeRefFrame!=NULL)
    activeReferenceOriginalValue = activeRefFrame->X;
  DFL=DFLin;
  PFL=PFLin;
  debug=debugIn;
  // Initialize action
  err = this->init();
  // Check for unhandled keywords
  actionArgs.CheckForMoreArgs();
  // Store the value of useMass set by the actions init
  useMassOriginalValue = useMass;
  // Store the value of useImage set by actions init
  useImageOriginalValue = useImage;

  return ( err );
}

// Action::Setup()
/** Set up action for the current parm file. This is where any 
  * parm-dependent variables should be set such as atom masks etc. The
  * current parm memory address is set here but can also be modified
  * by the action, this allows e.g. stripping of the parm. Only copies
  * of the parm should be modified; a reference to the original parm is
  * always stored in CpptrajState and can be reset there with the 'unstrip'
  * command.
  * \param ParmAddress memory address of current parm; may be changed
  *        by the action.
  */
int Action::Setup(AmberParm **ParmAddress) {
  int err;
  
  currentParm = *ParmAddress;
  // If useImage, check imaging type based on prmtop box.
  useImage = useImageOriginalValue;
  if (!useImage)
    imageType = NOBOX;
  else {
    imageType = currentParm->boxType;
    if (currentParm->boxType==NOBOX && debug>0) 
      mprintf("    Warning: No box info in %s, disabling imaging.\n",currentParm->parmName);
  }
  // If useMass, check that parm actually has masses.
  useMass = useMassOriginalValue;
  if (currentParm->mass==NULL && useMass) {
    mprintf("    Warning: %s: Mass for this parm is NULL.\n",actionArgs.Command());
    mprintf("             Geometric center will be used instead of center of mass.\n");
    useMass=false;
  }
  // If an active reference frame has been defined, check that the number of
  // atoms in the reference frame matches the number of atoms in the 
  // current parm. If less, distance based mask parsing could cause a segfault,
  // so dont allow. If more, print a warning that the reference coords may not
  // match up with current parm.
  activeReference = activeReferenceOriginalValue;
  if (activeReference!=NULL) {
    if (activeRefFrame->natom < currentParm->natom) {
      mprintf("Warning: # atoms in active reference frame (%i) < # atoms in current\n",
              activeRefFrame->natom);
      mprintf("         topology %s (%i). Distance-based mask parsing will not work.\n",
              currentParm->parmName, currentParm->natom);
      activeReference = NULL;
    } else if (activeRefFrame->natom > currentParm->natom) {
      mprintf("Warning: # atoms in active reference frame (%i) > # atoms in current\n",
              activeRefFrame->natom);
      mprintf("         topology %s (%i). Distance-based mask parsing may not work.\n",
              currentParm->parmName, currentParm->natom);
      mprintf("         Ensure that there is a 1 to 1 correspondance between active\n");
      mprintf("         reference frame and topology %s\n",currentParm->parmName);
    }
  }
  // Set up actions for this parm
  err = this->setup();
  if (err) return err;
  // Set the value of parm address in case parm was changed, e.g. in strip
  *ParmAddress = currentParm;
  return 0;
}

// Action::DoAction() 
/** Perform action on the current frame. The current frame memory address
  * is passed in and can be modified by the action, again for things like
  * stripping etc. The current frame number is also passed in.
  * \param FrameAddress memory address of the current frame; may be
  *        changed by the action.
  * \param frameNumIn number of the current frame
  */
Action::ActionReturnType Action::DoAction(Frame **FrameAddress, int frameNumIn) {
  ActionReturnType err;

  currentFrame = *FrameAddress;
  frameNum = frameNumIn;
  err = (ActionReturnType)this->action(); // NOTE: Fix return type eventually
  // Any state but ok means do not modify the frame. Return now.
  if (err!=ACTION_OK) return err;
  // Set the value of frame address in case frame was changed, e.g. in strip
  *FrameAddress = currentFrame;
  return ACTION_OK;
}


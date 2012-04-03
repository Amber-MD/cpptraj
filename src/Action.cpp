#include "Action.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action::Action() :
  noInit(false), 
  noSetup(false), 
  isSeparate(false),
  currentParm(NULL),
  currentFrame(NULL),
  DSL(NULL),
  DFL(NULL),
  PFL(NULL),
  FL(NULL),
  useMass(false),
  useMassOriginalValue(false),
  useImage(false),
  useImageOriginalValue(false),
  imageType (NOIMAGE),
  debug(0),
  frameNum(0) 
{}

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
  * \param PFLin pointer to the master TopologyList
  * \param debugIn Debug level that action should be set to
  */
int Action::Init(DataSetList *DSLin, FrameList *FLin, DataFileList *DFLin, 
                 TopologyList *PFLin, int debugIn) 
{
  DSL=DSLin;
  FL=FLin;
  DFL=DFLin;
  PFL=PFLin;
  debug=debugIn;
  // Initialize action
  int err = this->init();
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
  * always stored in Cpptraj and can be reset there with the 'unstrip'
  * command.
  * \param ParmAddress memory address of current parm; may be changed
  *        by the action.
  */
int Action::Setup(Topology **ParmAddress) {
  currentParm = *ParmAddress;
  // Set imaging to value set by init() 
  useImage = useImageOriginalValue;
  if (!useImage)
    // Imaging disabled
    imageType = NOIMAGE;
  else {
    // Set imaging based on parm box.
    Box::BoxType parmboxtype = currentParm->BoxType();
    if (parmboxtype == Box::NOBOX) {
      imageType = NOIMAGE;
      if (debug>0)
        mprintf("    Warning: No box info in %s, disabling imaging.\n",currentParm->c_str());
    } else if (parmboxtype == Box::ORTHO)
      imageType = ORTHO;
    else 
      imageType = NONORTHO;
  }
  // If useMass, check that parm actually has masses.
  // NOTE: Mass is now always set to 1 if not read in so this only depends
  //       on what the action set useMass to.
  useMass = useMassOriginalValue;
/*  if (currentParm->mass==NULL && useMass) {
    mprintf("    Warning: %s: Mass for this parm is NULL.\n",actionArgs.Command());
    mprintf("             Geometric center will be used instead of center of mass.\n");
    useMass=false;
  }*/
  // Set up actions for this parm
  int err = this->setup();
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
Action::ActionReturnType Action::DoAction(Frame **FrameAddress, int frameNumIn) 
{
  currentFrame = *FrameAddress;
  frameNum = frameNumIn;
  ActionReturnType err = (ActionReturnType)this->action(); // NOTE: Fix return type eventually
  // Any state but ok means do not modify the frame. Return now.
  if (err!=ACTION_OK) return err;
  // Set the value of frame address in case frame was changed, e.g. in strip
  *FrameAddress = currentFrame;
  return ACTION_OK;
}


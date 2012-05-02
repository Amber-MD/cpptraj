// ActionList
#include "ActionList.h"
#include "CpptrajStdio.h"
// All action classes go here
// TODO: Re-enable actions
#include "Action_Distance.h"
#include "Action_Rmsd.h"
#include "Action_Dihedral.h"
#include "Action_Angle.h"
#include "Action_AtomMap.h"
#include "Action_Strip.h"
#include "Action_DSSP.h"
#include "Action_Center.h"
#include "Action_Hbond.h"
#include "Action_Image.h"
#include "Action_Surf.h"
#include "Action_Radgyr.h"
#include "Action_Mask.h"
#include "Action_Closest.h"
#include "Action_NAstruct.h"
#include "Action_Pucker.h"
#include "Action_Outtraj.h"
#include "Action_Rms2d.h"
#include "Action_Average.h"
#include "Action_Radial.h"
#include "Action_DistRmsd.h"
#include "Action_Jcoupling.h"
#include "Action_Clustering.h"
//#include "Action_Pairwise.h"
#include "Action_PtrajAction.h"
#include "Action_Molsurf.h"
#include "Action_CheckStructure.h"
//#include "Action_DihedralScan.h"
#include "Action_Rotdif.h"
#include "Action_RunningAvg.h"
#include "Action_RmsAvgCorr.h"
#include "Action_AtomicFluct.h"
#include "Action_Watershell.h"
#include "Action_AvgCoord.h"
#include "Action_Contacts.h"

// CONSTRUCTOR
ActionList::ActionList() :
  debug_(0)
{}

// DESTRUCTOR
ActionList::~ActionList() {
    // No need to cast back to whatever action was allocd since Action destructor is virtual
    for (action_it act = actionlist_.begin(); act != actionlist_.end(); ++act)
      delete *act; 
}

// ActionList::SetDebug()
void ActionList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0)
    mprintf("ActionList DEBUG LEVEL SET TO %i\n",debug_);
}

// ActionList::AddAction()
/** Check if the first argument of the given arglist is an action keyword.
  * if so set up the appropriate action class.
  * \param argIn input argument list
  * \return 0 if action successfully added to the list, 1 if argument
  *           not recognized.
  */
int ActionList::AddAction(ArgList &argIn) {
  Action *Act;

  // Decide what action this is based on the command.
  if      (argIn.CommandIs("distance"))       {Act=new Distance;}
  else if (argIn.CommandIs("rms2d"))          {Act=new Rms2d;   }
  else if (argIn.CommandIs("2drms"))          {Act=new Rms2d;   }
  else if (argIn.CommandIs("rmsd"))           {Act=new Rmsd;    }
  else if (argIn.CommandIs("rms"))            {Act=new Rmsd;    }
  else if (argIn.CommandIs("dihedral"))       {Act=new Dihedral;}
  else if (argIn.CommandIs("atommap"))        {Act=new Action_AtomMap; }
  else if (argIn.CommandIs("angle"))          {Act=new Angle;   }
  else if (argIn.CommandIs("strip"))          {Act=new Strip;   }
  else if (argIn.CommandIs("secstruct"))      {Act=new DSSP;    }
  else if (argIn.CommandIs("center"))         {Act=new Center;  }
  else if (argIn.CommandIs("hbond"))          {Act=new Hbond;   }
  else if (argIn.CommandIs("image"))          {Act=new Image;   }
  else if (argIn.CommandIs("surf"))           {Act=new Surf;    }
  else if (argIn.CommandIs("radgyr"))         {Act=new Radgyr;  }
  else if (argIn.CommandIs("mask"))           {Act=new ActionMask;}
  else if (argIn.CommandIs("closest"))        {Act=new Closest; }
  else if (argIn.CommandIs("nastruct"))       {Act=new NAstruct;}
  else if (argIn.CommandIs("pucker"))         {Act=new Pucker;  }
  else if (argIn.CommandIs("outtraj"))        {Act=new Outtraj; }
  else if (argIn.CommandIs("unstrip"))        {Act=new Unstrip; }
  else if (argIn.CommandIs("average"))        {Act=new Average; }
  else if (argIn.CommandIs("radial"))         {Act=new Radial;  }
  else if (argIn.CommandIs("drmsd"))          {Act=new DistRmsd;}
  else if (argIn.CommandIs("drms"))           {Act=new DistRmsd;}
  else if (argIn.CommandIs("jcoupling"))      {Act=new Jcoupling;}
  else if (argIn.CommandIs("cluster"))        {Act=new Clustering;}
  //else if (argIn.CommandIs("pairwise"))       {Act=new Pairwise;}
  else if (argIn.CommandIs("molsurf"))        {Act=new Molsurf; }
  else if (argIn.CommandIs("checkstructure")) {Act=new CheckStructure;}
  else if (argIn.CommandIs("check"))          {Act=new CheckStructure;}
  //else if (argIn.CommandIs("dihedralscan"))   {Act=new DihedralScan;}
  else if (argIn.CommandIs("rotdif"))         {Act=new Rotdif;}
  else if (argIn.CommandIs("runningaverage")) {Act=new RunningAvg;}
  else if (argIn.CommandIs("runavg"))         {Act=new RunningAvg;}
  else if (argIn.CommandIs("runningaverage")) {Act=new RunningAvg;}
  else if (argIn.CommandIs("rmsavgcorr"))     {Act=new RmsAvgCorr;}
  else if (argIn.CommandIs("atomicfluct"))    {Act=new AtomicFluct;}
  else if (argIn.CommandIs("watershell"))     {Act=new Watershell;}
  else if (argIn.CommandIs("avgcoord"))       {Act=new Action_AvgCoord;}
  else if (argIn.CommandIs("contactstest"))   {Act=new Action_Contacts;}
  // PTRAJ
  else if (argIn.CommandIs("atomicfluct3D") ||
           argIn.CommandIs("contacts") ||
           argIn.CommandIs("correlation") ||
           argIn.CommandIs("clusterdihedral") ||
           argIn.CommandIs("diffusion") ||
           argIn.CommandIs("dipole") ||
           argIn.CommandIs("dnaiontracker") ||
           argIn.CommandIs("echo") ||
           argIn.CommandIs("grid") ||
           argIn.CommandIs("matrix") ||
           argIn.CommandIs("principal") ||
           argIn.CommandIs("projection") ||
           argIn.CommandIs("randomizeions") ||
           argIn.CommandIs("scale") ||
           argIn.CommandIs("unwrap") ||
           argIn.CommandIs("vector") 
          )
  {
    Act = new PtrajAction;
  } else return 1; 

  // Pass in the argument list
  Act->SetArg(argIn);
  // Debug
  if (debug_>0) mprintf("    Added action %s\n", Act->ActionCommand());

  // Store action in list
  actionlist_.push_back(Act);

  return 0;  
}

// ActionList::Init()
/** Initialize non-parm-specific data for each action (like datasets). If an 
  * action cannot be initialized deactivate it. Also set action debug level.
  */
int ActionList::Init( DataSetList *DSL, FrameList *FL, DataFileList *DFL, 
                      TopologyList *PFL, bool exitOnError) 
{
  mprintf("\nACTIONS: Initializing %zu actions:\n",actionlist_.size());
  unsigned int actnum = 0;
  for (action_it act = actionlist_.begin(); act != actionlist_.end(); ++act)
  {
    mprintf("  %u: [%s]\n", actnum++, (*act)->CmdLine());
    if ((*act)->NoInit()) {
      mprintf("Warning: Action %s is not active.\n", (*act)->ActionCommand());
    } else {
      if ( (*act)->Init( DSL, FL, DFL, PFL, debug_ ) ) {
        if (exitOnError) {
          mprinterr("Error: Init failed for [%s].\n", (*act)->CmdLine());
          return 1;
        } else {
          mprintf("Warning: Init failed for [%s]: DEACTIVATING\n",
                  (*act)->CmdLine());
          (*act)->SetNoInit();
        }
      }
    }
    mprintf("\n");
  }

  return 0;
}

// ActionList::Setup()
/** Attempt to set up all actions in the action list with the given parm
  * If an action cannot be set up skip it.
  */
int ActionList::Setup(Topology **ParmAddress) {
  Topology *OriginalParm = *ParmAddress;

  mprintf(".....................................................\n");
  mprintf("PARM [%s]: Setting up %zu actions.\n",(*ParmAddress)->c_str(),actionlist_.size());
  unsigned int actnum = 0;
  for (action_it act = actionlist_.begin(); act != actionlist_.end(); ++act)
  {
    if ((*act)->NoInit()) {
      // If action was not initialized, it cannot be setup.
      (*act)->SetNoSetup();
    } else {
      // Attempt to set up action.
      mprintf("  %u: [%s]\n", actnum++, (*act)->CmdLine());
      (*act)->ResetSetup();
      int err = (*act)->Setup(ParmAddress);
      if (err==1) {
        mprintf("Warning: Setup failed for [%s]: Skipping\n",
                (*act)->CmdLine());
        (*act)->SetNoSetup();
        //return 1;
      } else if (err==2) {
        // Return value of 2 requests return to original parm
        *ParmAddress = OriginalParm;
      }
      //fprintf(stdout,"DEBUG: After Action %i Setup parmName is %s\n",act,P->parmName);
    }
  }
  //mprintf(".....................................................\n");

  return 0;
}

// ActionList::DoActions()
/** Perform actions in the action list on the given Frame. Skip actions not 
  * initialized or not setup. 
  * \param FrameAddress Memory address of the current frame.
  * \param frameNumIn The current frame number.
  * \return true if coordinate output should be suppressed.
  * \return false if coordinate output should be performed.
  */
bool ActionList::DoActions(Frame **FrameAddress, int frameNumIn) {
  Frame *OriginalFrame = *FrameAddress;

  //fprintf(stdout,"DEBUG: Performing %i actions on frame %i.\n",Naction,frameNumIn);
  for (action_it act = actionlist_.begin(); act != actionlist_.end(); ++act) 
  {
    // Skip deactivated actions
    if ((*act)->NoSetup()) continue;
    // Perform action on frame
    Action::ActionReturnType err = (*act)->DoAction(FrameAddress, frameNumIn);
    // Check for action special conditions/errors
    if (err != Action::ACTION_OK) {
      if (err == Action::ACTION_USEORIGINALFRAME) {
        // Return value of 2 requests return to original frame
        *FrameAddress = OriginalFrame;
      } else if (err == Action::ACTION_SUPPRESSCOORDOUTPUT) {
        // Skip the rest of the actions and suppress output. Necessary when
        // e.g. performing a running average over coords.
        return true;
      } else {
        // If here return type is ACTION_ERR.
        // Treat actions that fail as if they could not be set up
        mprintf("Warning: Action [%s] failed, frame %i.\n", (*act)->CmdLine(),
              frameNumIn);
        (*act)->SetNoSetup();
      }
    } 
  }
  return false;
}

// ActionList::Print()
void ActionList::Print() {
  mprintf("\nACTION OUTPUT:\n");
  for (action_it act = actionlist_.begin(); act != actionlist_.end(); ++act)
  {
    // Skip deactivated actions
    if ((*act)->NoInit()) continue;
    (*act)->print();
  }
}

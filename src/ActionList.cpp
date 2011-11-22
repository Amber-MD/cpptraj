// ActionList
#include "ActionList.h"
#include "CpptrajStdio.h"
// All action classes go here
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
#include "Action_Pairwise.h"
#include "Action_PtrajAction.h"
#include "Action_Molsurf.h"
#include "Action_CheckStructure.h"
#include "Action_DihedralScan.h"

// CONSTRUCTOR
ActionList::ActionList() {
  Naction=0;
  debug=0;
}

// DESTRUCTOR
ActionList::~ActionList() {
    // No need to cast back to whatever action was allocd since Action destructor is virtual
    for (int i=0; i<Naction; i++) 
      delete actionlist[i];
}

// ActionList::SetDebug()
void ActionList::SetDebug(int debugIn) {
  debug=debugIn;
  if (debug>0)
    mprintf("ActionList DEBUG LEVEL SET TO %i\n",debug);
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
  else if (argIn.CommandIs("atommap"))        {Act=new AtomMap; }
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
  else if (argIn.CommandIs("pairwise"))       {Act=new Pairwise;}
  else if (argIn.CommandIs("molsurf"))        {Act=new Molsurf; }
  else if (argIn.CommandIs("checkstructure")) {Act=new CheckStructure;}
  else if (argIn.CommandIs("check"))          {Act=new CheckStructure;}
  else if (argIn.CommandIs("dihedralscan"))   {Act=new DihedralScan;}
  // PTRAJ
  else if (argIn.CommandIs("atomicfluct") ||
           argIn.CommandIs("atomicfluct3D") ||
           argIn.CommandIs("checkoverlap") ||
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
           argIn.CommandIs("runningaverage") ||
           argIn.CommandIs("scale") ||
           argIn.CommandIs("unwrap") ||
           argIn.CommandIs("vector") ||
           argIn.CommandIs("watershell") )
  {
    Act = new PtrajAction;
  } else return 1; 

  // Pass in the argument list
  Act->SetArg(argIn);
  // Debug
  if (debug>0) mprintf("    Added action %s\n", Act->ActionCommand());

  // Store action in list
  actionlist.push_back(Act);
  Naction++;

  return 0;  
}

// ActionList::Init()
/** Initialize non-parm-specific data for each action (like datasets). If an 
  * action cannot be initialized deactivate it. Also set action debug level.
  */
int ActionList::Init( DataSetList *DSL, FrameList *FL, 
                           DataFileList *DFL, ParmFileList *PFL) {
  mprintf("\nACTIONS: Initializing %i actions:\n",Naction);
  for (int act=0; act<Naction; act++) {
    mprintf("  %i: [%s]\n",act,actionlist[act]->CmdLine());
    if (actionlist[act]->noInit) {
      mprintf("    WARNING: Action %s is not active.\n",actionlist[act]->ActionCommand());
    } else {
      if ( actionlist[act]->Init( DSL, FL, DFL, PFL, debug ) ) {
        mprintf("    WARNING: Init failed for [%s]: DEACTIVATING\n",
                actionlist[act]->CmdLine());
        actionlist[act]->noInit=true;
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
int ActionList::Setup(AmberParm **ParmAddress) {
  int err;
  AmberParm *OriginalParm = *ParmAddress;

  mprintf("  ...... Setting up %i actions for %s ....\n",Naction,(*ParmAddress)->parmName);
  for (int act=0; act<Naction; act++) {
    if (!actionlist[act]->noInit) {
      mprintf("  %i: [%s]\n",act,actionlist[act]->CmdLine());
      actionlist[act]->noSetup=false;
      err = actionlist[act]->Setup(ParmAddress);
      if (err==1) {
        mprintf("      WARNING: Setup failed for [%s]: Skipping\n",
                actionlist[act]->CmdLine());
        actionlist[act]->noSetup=true;
        //return 1;
      } else if (err==2) {
        // Return value of 2 requests return to original parm
        *ParmAddress = OriginalParm;
      }
      //fprintf(stdout,"DEBUG: After Action %i Setup parmName is %s\n",act,P->parmName);
    }
  }
  mprintf("  .....................................................\n");

  return 0;
}

// ActionList::DoActions()
/** Perform actions in the action list on the given Frame. Skip actions not 
  * initialized or not setup. frameNumIn is the current frame number.
  */
void ActionList::DoActions(Frame **FrameAddress, int frameNumIn) {
  int err;
  Frame *OriginalFrame = *FrameAddress;

  //fprintf(stdout,"DEBUG: Performing %i actions on frame %i.\n",Naction,frameNumIn);
  for (int act=0; act<Naction; act++) {
    // Skip deactivated actions
    if (actionlist[act]->noInit || actionlist[act]->noSetup) continue;
    err = actionlist[act]->DoAction(FrameAddress, frameNumIn);
    if (err==1) {
      // Treat actions that fail as if they could not be set up
      mprintf("Warning: Action [%s] failed, frame %i.\n",actionlist[act]->CmdLine(),
              frameNumIn);
      actionlist[act]->noSetup=true;
    } else if (err==2) {
      // Return value of 2 requests return to original frame
      *FrameAddress = OriginalFrame;
    }
  }
}

// ActionList::Print()
void ActionList::Print() {
  mprintf("\nACTION OUTPUT:\n");
  for (int act=0; act<Naction; act++) {
    // Skip deactivated actions
    if (actionlist[act]->noInit) continue;
    actionlist[act]->print();
  }
}

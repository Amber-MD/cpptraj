// ActionList
#include <cstdlib>
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

// CONSTRUCTOR
ActionList::ActionList() {
  actionlist=NULL;
  Naction=0;
  debug=0;
}

// DESTRUCTOR
ActionList::~ActionList() {
  if (actionlist!=NULL) {
    // No need to cast back to whatever action was allocd since Action destructor is virtual
    for (int i=0; i<Naction; i++) 
      delete actionlist[i];
    free(actionlist);
  } 
}

/* ActionList::SetDebug()
 * Set Action debug level
 */
void ActionList::SetDebug(int debugIn) {
  debug=debugIn;
  if (debug>0)
    mprintf("ActionList DEBUG LEVEL SET TO %i\n",debug);
}

/* ActionList::Add()
 * Add a specific type of action class to the action list. 
 */
int ActionList::Add(ArgList *A) {
  Action *Act;

  // Decide what action this is based on the command.
  if      (A->CommandIs("distance")) {Act=new Distance;}
  else if (A->CommandIs("rms2d"))    {Act=new Rms2d;   }
  else if (A->CommandIs("2drms"))    {Act=new Rms2d;   }
  else if (A->CommandIs("rmsd",3))   {Act=new Rmsd;    }
  else if (A->CommandIs("dihedral")) {Act=new Dihedral;}
  else if (A->CommandIs("atommap"))  {Act=new AtomMap; }
  else if (A->CommandIs("angle"))    {Act=new Angle;   }
  else if (A->CommandIs("strip"))    {Act=new Strip;   }
  else if (A->CommandIs("secstruct")){Act=new DSSP;    }
  else if (A->CommandIs("center"))   {Act=new Center;  }
  else if (A->CommandIs("hbond"))    {Act=new Hbond;   }
  else if (A->CommandIs("image"))    {Act=new Image;   }
  else if (A->CommandIs("surf"))     {Act=new Surf;    }
  else if (A->CommandIs("radgyr"))   {Act=new Radgyr;  }
  else if (A->CommandIs("mask"))     {Act=new ActionMask;}
  else if (A->CommandIs("closest"))  {Act=new Closest; }
  else if (A->CommandIs("nastruct")) {Act=new NAstruct;}
  else if (A->CommandIs("pucker"))   {Act=new Pucker;  }
  else if (A->CommandIs("outtraj"))  {Act=new Outtraj; }
  else if (A->CommandIs("unstrip"))  {Act=new Unstrip; }
  else if (A->CommandIs("average"))  {Act=new Average; }
  else if (A->CommandIs("radial"))   {Act=new Radial;  }
  else if (A->CommandIs("drmsd",4))  {Act=new DistRmsd;}
  else if (A->CommandIs("jcoupling")){Act=new Jcoupling;}
  else if (A->CommandIs("cluster"))  {Act=new Clustering;}
  else return 1; 

  // Pass in the argument list
  Act->setArg(A);
  // Debug
  if (debug>0) mprintf("    Added action %s\n", Act->Name());

  // Store action in list
  actionlist=(Action**) realloc(actionlist,(Naction+1) * sizeof(Action*));
  actionlist[Naction]=Act;
  Naction++;

  return 0;  
}

/* ActionList::Init()
 * Initialize non-parm-specific data for each action (like datasets). If an 
 * action cannot be initialized deactivate it. Also set action debug level.
 */
int ActionList::Init( DataSetList *DSL, FrameList *FL, 
                           DataFileList *DFL, ParmFileList *PFL) {
  mprintf("\nACTIONS: Initializing %i actions:\n",Naction);
  for (int act=0; act<Naction; act++) {
    mprintf("  %i: [%s]\n",act,actionlist[act]->CmdLine());
    if (actionlist[act]->noInit) {
      mprintf("    WARNING: Action %s is not active.\n",actionlist[act]->Name());
    } else {
      actionlist[act]->ResetArg();
      if ( actionlist[act]->Init( DSL, FL, DFL, PFL, debug ) ) {
        mprintf("    WARNING: Init failed for [%s]: DEACTIVATING\n",
                actionlist[act]->CmdLine());
        actionlist[act]->noInit=1;
      }
    }
    mprintf("\n");
  }

  return 0;
}

/* ActionList::Setup()
 * Attempt to set up all actions in the action list with the given parm
 * If an action cannot be set up skip it.
 */
int ActionList::Setup(AmberParm **ParmAddress) {
  int err;
  AmberParm *OriginalParm = *ParmAddress;

  mprintf("    .... Setting up %i actions for %s ....\n",Naction,(*ParmAddress)->parmName);
  for (int act=0; act<Naction; act++) {
    if (actionlist[act]->noInit==0) {
      if (debug>0) mprintf("    %i: [%s]\n",act,actionlist[act]->CmdLine());
      actionlist[act]->noSetup=0;
      actionlist[act]->ResetArg();
      err = actionlist[act]->Setup(ParmAddress);
      if (err==1) {
        mprintf("      WARNING: Setup failed for [%s]: Skipping\n",
                actionlist[act]->CmdLine());
        actionlist[act]->noSetup=1;
        //return 1;
      } else if (err==2) {
        // Return value of 2 requests return to original parm
        *ParmAddress = OriginalParm;
      }
      //fprintf(stdout,"DEBUG: After Action %i Setup parmName is %s\n",act,P->parmName);
    }
  }
  mprintf("    ...................................................\n");

  return 0;
}

/* ActionList::DoActions()
 * Perform actions in the action list on the given Frame. Skip actions not 
 * initialized or not setup. frameIn is the current frame number.
 */
void ActionList::DoActions(Frame **FrameAddress, int frameIn) {
  int err;
  Frame *OriginalFrame = *FrameAddress;

  //fprintf(stdout,"DEBUG: Performing %i actions on frame %i.\n",Naction,frameIn);
  for (int act=0; act<Naction; act++) {
    // Skip deactivated actions
    if (actionlist[act]->noInit || actionlist[act]->noSetup) continue;
    err = actionlist[act]->DoAction(FrameAddress, frameIn);
    if (err==1) {
      // Treat actions that fail as if they could not be set up
      mprintf("Warning: Action [%s] failed, frame %i.\n",actionlist[act]->CmdLine(),
              frameIn);
      actionlist[act]->noSetup=1;
    } else if (err==2) {
      // Return value of 2 requests return to original frame
      *FrameAddress = OriginalFrame;
    }
  }
}

/* ActionList::Print()
 * Non-dataset print for all actions. 
 */
void ActionList::Print() {
  mprintf("\nACTION OUTPUT:\n");
  for (int act=0; act<Naction; act++) {
    // Skip deactivated actions
    if (actionlist[act]->noInit) continue;
    actionlist[act]->print();
  }
}

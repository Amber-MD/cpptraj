// PtrajActionList
#include <cstdlib>
#include "PtrajActionList.h"
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

// Constructor
PtrajActionList::PtrajActionList() {
  ActionList=NULL;
  Naction=0;
  debug=0;
}

// Destructor
PtrajActionList::~PtrajActionList() {
  int i;

  if (ActionList!=NULL) {
    // No need to cast back to whatever action was allocd since Action destructor is virtual
    for (i=0; i<Naction; i++) 
      delete ActionList[i];
    free(ActionList);
  } 
}

// Set Action debug level
void PtrajActionList::SetDebug(int debugIn) {
  debug=debugIn;
  if (debug>0)
    mprintf("PtrajActionList DEBUG LEVEL SET TO %i\n",debug);
}

/*
 * PtrajActionList::Add()
 * Add a specific type of action class to the action list. 
 */
int PtrajActionList::Add(ArgList *A) {
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
  else return 1; 

  // Pass in the argument list
  Act->setArg(A);
  // Debug
  if (debug>0) mprintf("    Added action %s\n", Act->Name());

  // Store action in list
  ActionList=(Action**) realloc(ActionList,(Naction+1) * sizeof(Action*));
  ActionList[Naction]=Act;
  Naction++;

  return 0;  
}

/* 
 * PtrajActionList::Init()
 * Initialize non-parm-specific data for each action (like datasets). If an 
 * action cannot be initialized deactivate it. Also set action debug level.
 */
int PtrajActionList::Init( DataSetList *DSL, FrameList *FL, 
                           DataFileList *DFL, ParmFileList *PFL) {
 int act;

  mprintf("\nACTIONS: Initializing %i actions:\n",Naction);
  for (act=0; act<Naction; act++) {
    mprintf("  %i: [%s]\n",act,ActionList[act]->CmdLine());
    if (ActionList[act]->noInit) {
      mprintf("    WARNING: Action %s is not active.\n",ActionList[act]->Name());
    } else {
      ActionList[act]->ResetArg();
      if ( ActionList[act]->Init( DSL, FL, DFL, PFL, debug ) ) {
        mprintf("    WARNING: Init failed for [%s]: DEACTIVATING\n",
                ActionList[act]->CmdLine());
        ActionList[act]->noInit=1;
      }
    }
    mprintf("\n");
  }

  return 0;
}


/* 
 * PtrajActionList::Setup()
 * Attempt to set up all actions in the action list with the given parm
 * If an action cannot be set up skip it.
 */
int PtrajActionList::Setup(AmberParm **ParmAddress) {
  int act,err;
  AmberParm *OriginalParm = *ParmAddress;

  mprintf("    .... Setting up %i actions for %s ....\n",Naction,(*ParmAddress)->parmName);
  for (act=0; act<Naction; act++) {
    if (ActionList[act]->noInit==0) {
      if (debug>0) mprintf("    %i: [%s]\n",act,ActionList[act]->CmdLine());
      ActionList[act]->noSetup=0;
      ActionList[act]->ResetArg();
      err = ActionList[act]->Setup(ParmAddress);
      if (err==1) {
        mprintf("      WARNING: Setup failed for [%s]: Skipping\n",
                ActionList[act]->CmdLine());
        ActionList[act]->noSetup=1;
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

/* 
 * PtrajActionList::DoActions()
 * Perform actions in the action list on the given Frame. Skip actions not 
 * initialized or not setup. frameIn is the current frame number.
 */
void PtrajActionList::DoActions(Frame **FrameAddress, int frameIn) {
  int act,err;
  Frame *OriginalFrame = *FrameAddress;

  //fprintf(stdout,"DEBUG: Performing %i actions on frame %i.\n",Naction,frameIn);
  for (act=0; act<Naction; act++) {
    // Skip deactivated actions
    if (ActionList[act]->noInit || ActionList[act]->noSetup) continue;
    err = ActionList[act]->DoAction(FrameAddress, frameIn);
    if (err==1) {
      // Treat actions that fail as if they could not be set up
      mprintf("Warning: Action [%s] failed, frame %i.\n",ActionList[act]->CmdLine(),
              frameIn);
      ActionList[act]->noSetup=1;
    } else if (err==2) {
      // Return value of 2 requests return to original frame
      *FrameAddress = OriginalFrame;
    }
  }

}

// Non-dataset print for all actions. 
void PtrajActionList::Print() {
  int act;

  mprintf("\nOUTPUT:\n");
  for (act=0; act<Naction; act++) {
    // Skip deactivated actions
    if (ActionList[act]->noInit) continue;
    ActionList[act]->print();
  }
}

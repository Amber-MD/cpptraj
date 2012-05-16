// PtrajAction
#include <cstdlib>
#include <cstring> // memset
#include "Action_PtrajAction.h"
#include "CpptrajStdio.h"
#include "ptraj_convert.h"

// CONSTRUCTOR
PtrajAction::PtrajAction() {
  //fprintf(stderr,"PtrajAction Con\n");
  actioninfo = NULL;
  argumentStack = NULL;
  x_coord = NULL;
  y_coord = NULL;
  z_coord = NULL;
  CalledSetup = false;
  ptraj_box[0] = 0;
  ptraj_box[1] = 0;
  ptraj_box[2] = 0;
  ptraj_box[3] = 0;
  ptraj_box[4] = 0;
  ptraj_box[5] = 0;
  coordinate_update = false;
} 

// DESTRUCTOR
PtrajAction::~PtrajAction() {
  //fprintf(stderr,"PtrajAction Destructor.\n");
  // Free reference information
  FreeReferenceInfo();
  // Free arguments
  FreeArgumentStack(argumentStack);
  // Free actionInfo
  if (actioninfo != NULL) {
    // Call actions cleanup routine
    if (CalledSetup)
      actioninfo->fxn(actioninfo, x_coord, y_coord, z_coord, ptraj_box, PTRAJ_CLEANUP);
    // Free state
    FreeState(actioninfo->state);
    free(actioninfo);
  }
  // Free coords
  if (x_coord!=NULL) free(x_coord);
  if (y_coord!=NULL) free(y_coord);
  if (z_coord!=NULL) free(z_coord);
}

// PtrajAction::init()
/** Set up the action information structure based on the given Ptraj 
  * command. Set the action type and function pointer to the corresponding
  * routine in ptraj_actions.c. Convert the argument list to a format
  * that is readable by the routines in ptraj_actions.c. Initialize
  * state memory.
  */
int PtrajAction::init( ) {
  // Set up action information structure
  actioninfo = (actionInformation *) malloc(sizeof(actionInformation));
  // DEBUG
  //mprintf("DEBUG:\taction address (init): %x\n",actioninfo);
  INITIALIZE_actionInformation(actioninfo);
  actioninfo->type = TRANSFORM_NOOP;

  // Set the action type and function based on the command
  /*if ( actionArgs.CommandIs("correlation")    ) {
    actioninfo->type = TRANSFORM_CORRELATION;
    actioninfo->fxn  = (actionFunction) transformCorr;
  } else*/ if ( actionArgs.CommandIs("clusterdihedral")) {
    actioninfo->type = TRANSFORM_DIHEDRALCLUSTER;
    actioninfo->fxn  = (actionFunction) transformDihedralCluster;
  } else if ( actionArgs.CommandIs("diffusion")      ) {
    actioninfo->type = TRANSFORM_DIFFUSION;
    actioninfo->fxn  = (actionFunction) transformDiffusion;
  } else if ( actionArgs.CommandIs("dipole")         ) {
    actioninfo->type = TRANSFORM_DIPOLE;
    actioninfo->fxn  = (actionFunction) transformDipole;
  } else if ( actionArgs.CommandIs("dnaiontracker")  ) {
    actioninfo->type = TRANSFORM_DNAIONTRACKER;
    actioninfo->fxn  = (actionFunction) transformDNAiontracker;
  } else if ( actionArgs.CommandIs("echo")           ) {
    actioninfo->type = TRANSFORM_ECHO;
    actioninfo->fxn = (actionFunction) transformEcho;
  } else if ( actionArgs.CommandIs("grid")           ) {
    actioninfo->type = TRANSFORM_GRID;
    actioninfo->fxn  = (actionFunction) transformGrid;
  } else if ( actionArgs.CommandIs("FEoH")           ) {
    actioninfo->type = TRANSFORM_FEOH;
    actioninfo->fxn  = (actionFunction) transformGibbsEnergyOfHydration;
  } else if ( actionArgs.CommandIs("principal")      ) {
    actioninfo->type = TRANSFORM_PRINCIPAL;
    actioninfo->fxn  = (actionFunction) transformPrincipal;
    coordinate_update = true;
  } else if ( actionArgs.CommandIs("projection")     ) {
    actioninfo->type = TRANSFORM_PROJECTION;
    actioninfo->fxn  = (actionFunction) transformProjection;
  } else if ( actionArgs.CommandIs("randomizeions")  ) {
    actioninfo->type = TRANSFORM_RANDOMIZEIONS;
    actioninfo->fxn  = (actionFunction) transformRandomizeIons;
    coordinate_update = true;
  } else if ( actionArgs.CommandIs("scale")          ) {
    actioninfo->type = TRANSFORM_SCALE;
    actioninfo->fxn  = (actionFunction) transformScale;
  } else if ( actionArgs.CommandIs("unwrap")         ) {
    actioninfo->type = TRANSFORM_UNWRAP;
    actioninfo->fxn  = (actionFunction) transformUnwrap;
    coordinate_update = true;
  } else {
    mprinterr("Error: PtrajAction: Unrecognized Ptraj command: %s\n",actionArgs.Command());
    return 1;
  }

  // Convert the ArgList to the argStackType used by functions in ptraj_actions.c
  argumentStack = CreateArgumentStack(actionArgs,debug);
  actioninfo->carg1 = (void *) &argumentStack;
  //mprintf("DEBUG:\taction->carg1 address (init): %x\n",actioninfo->carg1);
  // NOTE: Should ptraj be freeing up the args?
  // DEBUG
  //argStackType **argumentStackPointer = (argStackType **) actioninfo->carg1;
  //mprintf("DEBUG:\targumentStack address (init 2): %x\n",*argumentStackPointer);

  // Set reference structure
  Frame *refframe = FL->GetFrame(0);
  if (refframe!=NULL) 
    SetReferenceInfo(refframe->CoordPtr(), refframe->Natom());

  // Dont call setup here since there is no state information yet
  mprintf("    PTRAJ ACTION: [%s]\n",actionArgs.Command());

  // Set prnlev
  SetPrnlev(debug);
  
  return 0;
}

// PtrajAction::setup()
/** Initialize the state for the action information structure. Call the
  * action routine with mode PTRAJ_SETUP.
  */
int PtrajAction::setup() {
  // DEBUG
  //mprintf("DEBUG:\taction address (setup): %x\n",actioninfo);
  //mprintf("DEBUG:\taction->carg1 address (setup): %x\n",actioninfo->carg1);
  //argStackType **argumentStackPointer = (argStackType **) actioninfo->carg1;
  //mprintf("DEBUG:\targumentStack address (setup): %x\n",*argumentStackPointer);
  //printArgumentStack(argumentStackPointer);
  actioninfo->state = CreateState(currentParm, DSL->MaxFrames());

  // Now that the state info has been allocated, call setup
  if (actioninfo->fxn(actioninfo, NULL, NULL, NULL, NULL, PTRAJ_SETUP) < 0) {
    mprinterr("Error: PtrajAction: Could not set up action %s\n",actionArgs.Command());
    return 1;
  }
  // Since the functions in ptraj_actions.c are not on the stack, they are
  // not automatically initialized, so dont want to call them in action, 
  // print, or cleanup if set up failed or wasnt called.
  CalledSetup = true;
  
  // Print status
  actioninfo->fxn(actioninfo, NULL, NULL, NULL, NULL, PTRAJ_STATUS);

  // Allocate space for the X Y and Z vectors
  if (x_coord!=NULL) free(x_coord);
  if (y_coord!=NULL) free(y_coord);
  if (z_coord!=NULL) free(z_coord);
  x_coord = (double*) malloc( currentParm->Natom() * sizeof(double));
  y_coord = (double*) malloc( currentParm->Natom() * sizeof(double));
  z_coord = (double*) malloc( currentParm->Natom() * sizeof(double));
        
  return 0;  
}

// PtrajAction::action()
int PtrajAction::action() {
  // NOTE: Not checking for Called Setup here since the noSetup
  // flag should catch it.
  // Convert coordinates array X0Y0Z0X1... into separate X Y Z vectors
  // Also protect state box coords
  currentFrame->ConvertToPtrajXYZ(x_coord, y_coord, z_coord, ptraj_box);

  // Ptraj actions return 1 on success, 0 on failure 
  if (actioninfo->fxn(actioninfo, x_coord, y_coord, z_coord, ptraj_box, PTRAJ_ACTION)==0)
    return 1;

  // If necessary, translate ptraj X Y Z vectors back to frame
  if (coordinate_update) 
    currentFrame->SetFromPtrajXYZ(x_coord, y_coord, z_coord);

  return 0;
} 

// PtrajAction::print()
void PtrajAction::print() {
  // NOTE: Is it ok to just call this will NULLs?
  if (CalledSetup)
    actioninfo->fxn(actioninfo, x_coord, y_coord, z_coord, ptraj_box, PTRAJ_PRINT);
}

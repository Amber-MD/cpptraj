// PtrajAction
#include <cstdlib>
#include <cstring>
#include "Action_PtrajAction.h"
#include "CpptrajStdio.h"
#include "ptraj_actions.h"

// CONSTRUCTOR
PtrajAction::PtrajAction() {
  //fprintf(stderr,"PtrajAction Con\n");
  actionptr = NULL;
  x_coord = NULL;
  y_coord = NULL;
  z_coord = NULL;
} 

// DESTRUCTOR
PtrajAction::~PtrajAction() {
  //fprintf(stderr,"PtrajAction Destructor.\n");
  if (actionptr != NULL) {
    actionInformation *tmpaction = (actionInformation*) actionptr;
    // Call actions cleanup routine
    tmpaction->fxn(tmpaction, x_coord, y_coord, z_coord, NULL, PTRAJ_CLEANUP);
    // Free args
    argStackType **argumentStackAddress = (argStackType **) tmpaction->carg1;
    argStackType *argumentStack = *argumentStackAddress;
    // free individual args?
    free(argumentStack->marked);
    free(argumentStack);
    // Free state
    free(tmpaction->state->solventMask);
    free(tmpaction->state->solventMoleculeStart);
    free(tmpaction->state->solventMoleculeStop);
    free(tmpaction->state);
    free(tmpaction);
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
  actionInformation *action = NULL;

  // Set up action information structure
  action = (actionInformation *) malloc(sizeof(actionInformation));
  INITIALIZE_actionInformation(action);
  action->type = TRANSFORM_NOOP;

  // Set the action type and function based on the command
  if      ( actionArgs.CommandIs("atomicfluct")    ) {
    action->type = TRANSFORM_ATOMICFLUCT;
    action->fxn = (actionFunction) transformAtomicFluct;
  } else if ( actionArgs.CommandIs("atomicfluct3D")  ) {
    action->type = TRANSFORM_ATOMICFLUCT3D;
    action->fxn = (actionFunction) transformAtomicFluct3D;
  } else if ( actionArgs.CommandIs("checkoverlap")   ) {
    action->type = TRANSFORM_CHECKOVERLAP;
    action->fxn  = (actionFunction) transformCheckOverlap;
  } else if ( actionArgs.CommandIs("contacts")       ) {
    action->type = TRANSFORM_CONTACTS;
    action->fxn  = (actionFunction) transformContacts;
  } else if ( actionArgs.CommandIs("correlation")    ) {
    action->type = TRANSFORM_CORRELATION;
    action->fxn  = (actionFunction) transformCorr;
  } else if ( actionArgs.CommandIs("clusterdihedral")) {
    action->type = TRANSFORM_DIHEDRALCLUSTER;
    action->fxn  = (actionFunction) transformDihedralCluster;
  } else if ( actionArgs.CommandIs("diffusion")      ) {
    action->type = TRANSFORM_DIFFUSION;
    action->fxn  = (actionFunction) transformDiffusion;
  } else if ( actionArgs.CommandIs("dipole")         ) {
    action->type = TRANSFORM_DIPOLE;
    action->fxn  = (actionFunction) transformDipole;
  } else if ( actionArgs.CommandIs("dnaiontracker")  ) {
    action->type = TRANSFORM_DNAIONTRACKER;
    action->fxn  = (actionFunction) transformDNAiontracker;
  } else if ( actionArgs.CommandIs("echo")           ) {
    action->type = TRANSFORM_ECHO;
    action->fxn = (actionFunction) transformEcho;
  } else if ( actionArgs.CommandIs("grid")           ) {
    action->type = TRANSFORM_GRID;
    action->fxn  = (actionFunction) transformGrid;
  } else if ( actionArgs.CommandIs("matrix")         ) {
    action->type = TRANSFORM_MATRIX;
    action->fxn  = (actionFunction) transformMatrix;
  } else if ( actionArgs.CommandIs("principal")      ) {
    action->type = TRANSFORM_PRINCIPAL;
    action->fxn  = (actionFunction) transformPrincipal;
  } else if ( actionArgs.CommandIs("projection")     ) {
    action->type = TRANSFORM_PROJECTION;
    action->fxn  = (actionFunction) transformProjection;
  } else if ( actionArgs.CommandIs("randomizeions")  ) {
    action->type = TRANSFORM_RANDOMIZEIONS;
    action->fxn  = (actionFunction) transformRandomizeIons;
  } else if ( actionArgs.CommandIs("runningaverage") ) {
    action->type = TRANSFORM_RUNNINGAVERAGE;
    action->fxn  = (actionFunction) transformRunningAverage;
  } else if ( actionArgs.CommandIs("scale")          ) {
    action->type = TRANSFORM_SCALE;
    action->fxn  = (actionFunction) transformScale;
  } else if ( actionArgs.CommandIs("unwrap")         ) {
    action->type = TRANSFORM_UNWRAP;
    action->fxn  = (actionFunction) transformUnwrap;
  } else if ( actionArgs.CommandIs("vector")         ) {
    action->type = TRANSFORM_VECTOR;
    action->fxn  = (actionFunction) transformVector;
  } else if ( actionArgs.CommandIs("watershell")     ) {
    action->type = TRANSFORM_WATERSHELL;
    action->fxn  = (actionFunction) transformWatershell;
  } else {
    mprinterr("Error: PtrajAction: Unrecognized Ptraj command: %s\n",actionArgs.Command());
    return 1;
  }

  // Convert the ArgList to the argStackType used by functions in ptraj_actions.c
  argStackType *argumentStack = (argStackType*) malloc( sizeof(argStackType) );
  int nargs = 0;
  char *currentArg = actionArgs.getNextString();
  while (currentArg != NULL) {
    argumentStack->arglist = (char**) realloc(argumentStack->arglist, (nargs+1) * sizeof(char*));
    argumentStack->arglist[nargs] = (char*) malloc( (strlen(currentArg)+1) * sizeof(char) );
    strcpy(argumentStack->arglist[nargs], currentArg);
    currentArg = actionArgs.getNextString();
    nargs++;
  }
  argumentStack->nargs = nargs;
  argumentStack->marked = (char*) malloc(nargs * sizeof(char));
  memset(argumentStack->marked, 'F', nargs);
  action->carg1 = (void *) &argumentStack;
  // NOTE: Should ptraj be freeing up the args?

  // Initialize state memory
  action->state = (ptrajState*) malloc( sizeof(ptrajState) );
  INITIALIZE_ptrajState( action->state );

  // Dont call setup here since there is no state information yet
  
  actionptr = (void*) action;
  return 0;
}

// PtrajAction::setup()
/** Initialize the state for the action information structure. Call the
  * action routine with mode PTRAJ_SETUP.
  */
int PtrajAction::setup() {
  actionInformation *action = (actionInformation*) actionptr;
  ptrajState *state = action->state;
  // Place a copy of the current state into the action
  state->box[0] = currentParm->Box[0];
  state->box[1] = currentParm->Box[1];
  state->box[2] = currentParm->Box[2];
  state->box[3] = currentParm->Box[3];
  state->box[4] = currentParm->Box[4];
  state->box[5] = currentParm->Box[5];
  state->masses = currentParm->mass;
  state->charges = currentParm->charge;
  state->atoms = currentParm->natom;
  state->residues = currentParm->nres;
  state->ipres = currentParm->resnums;
  state->IFBOX = AmberIfbox(currentParm->Box[4]);
  //state->boxfixed
  state->molecules = currentParm->molecules;
  state->moleculeInfo = currentParm->atomsPerMol;
  // Solvent info
  state->solventMask = (int*) malloc( currentParm->natom * sizeof(int));
  state->solventMolecules = currentParm->solventMolecules;
  state->solventMoleculeStart = (int*) malloc( currentParm->solventMolecules * sizeof(int));
  state->solventMoleculeStop = (int*) malloc( currentParm->solventMolecules * sizeof(int));
  // Convert solvent mask
  if (currentParm->solventMask!=NULL) {
    for (int atom = 0; atom < currentParm->natom; atom++) {
      if (currentParm->solventMask[atom]=='T')
        state->solventMask[atom] = 1;
      else
        state->solventMask[atom] = 0;
    }
    for (int mol = 0; mol < currentParm->solventMolecules; mol++) {
      state->solventMoleculeStart[mol] = currentParm->solventMoleculeStart[mol];
      state->solventMoleculeStop[mol] = currentParm->solventMoleculeStop[mol];
    }
  } else {
    memset(state->solventMask,0,currentParm->natom);
    memset(state->solventMoleculeStart,0,currentParm->solventMolecules);
    memset(state->solventMoleculeStop,0,currentParm->solventMolecules);
  }
  state->solventAtoms = state->solventMoleculeStop[currentParm->solventMolecules-1] - 
                        state->solventMoleculeStart[0];
  state->atomName = currentParm->names;
  state->residueName = currentParm->resnames;
  state->maxFrames = currentParm->parmFrames; // NOTE: Is this correct?
  state->temp0 = 0.0;

  // Now that the state info has been allocated, call setup
  if (action->fxn(action, NULL, NULL, NULL, NULL, PTRAJ_SETUP) < 0) {
    mprinterr("Error: PtrajAction: Could not set up action %s\n",actionArgs.Command());
    return 1;
  }

  // Print status
  action->fxn(action, NULL, NULL, NULL, NULL, PTRAJ_STATUS);

  // Allocate space for the X Y and Z vectors
  if (x_coord!=NULL) free(x_coord);
  if (y_coord!=NULL) free(y_coord);
  if (z_coord!=NULL) free(z_coord);
  x_coord = (double*) malloc( currentParm->natom * sizeof(double));
  y_coord = (double*) malloc( currentParm->natom * sizeof(double));
  z_coord = (double*) malloc( currentParm->natom * sizeof(double));
        
  return 0;  
}

// PtrajAction::action()
int PtrajAction::action() {
  actionInformation *action = (actionInformation*) actionptr;
  // Convert coordinates array X0Y0Z0X1... into separate X Y Z vectors
  int i3 = 0;
  for (int atom = 0; atom < currentFrame->natom; atom++) {
    x_coord[atom] = currentFrame->X[i3++];
    y_coord[atom] = currentFrame->X[i3++];
    z_coord[atom] = currentFrame->X[i3++];
  }
  
  action->fxn(action, x_coord, y_coord, z_coord, currentFrame->box, PTRAJ_ACTION);

  return 0;
} 

// PtrajAction::print()
void PtrajAction::print() {
  actionInformation *action = (actionInformation*) actionptr;
  // NOTE: Is it ok to just call this will NULLs?
  action->fxn(action, x_coord, y_coord, z_coord, currentParm->Box, PTRAJ_PRINT);
}

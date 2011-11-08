// PtrajAction
#include <cstdlib>
#include <cstring>
#include "Action_PtrajAction.h"
#include "CpptrajStdio.h"
#include "ptraj_actions.h"

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
} 

// DESTRUCTOR
PtrajAction::~PtrajAction() {
  //fprintf(stderr,"PtrajAction Destructor.\n");
  if (argumentStack!=NULL) {
    // free individual args?
    for (int arg=0; arg < argumentStack->nargs; arg++)
      free(argumentStack->arglist[arg]);
    free(argumentStack->arglist);
    free(argumentStack->marked);
    free(argumentStack);
  }
  if (actioninfo != NULL) {
    // Call actions cleanup routine
    if (CalledSetup)
      actioninfo->fxn(actioninfo, x_coord, y_coord, z_coord, ptraj_box, PTRAJ_CLEANUP);
    // Free state
    free(actioninfo->state->ipres);
    free(actioninfo->state->solventMask);
    free(actioninfo->state->solventMoleculeStart);
    free(actioninfo->state->solventMoleculeStop);
    free(actioninfo->state);
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
  if      ( actionArgs.CommandIs("atomicfluct")    ) {
    actioninfo->type = TRANSFORM_ATOMICFLUCT;
    actioninfo->fxn = (actionFunction) transformAtomicFluct;
  } else if ( actionArgs.CommandIs("atomicfluct3D")  ) {
    actioninfo->type = TRANSFORM_ATOMICFLUCT3D;
    actioninfo->fxn = (actionFunction) transformAtomicFluct3D;
  } else if ( actionArgs.CommandIs("checkoverlap")   ) {
    actioninfo->type = TRANSFORM_CHECKOVERLAP;
    actioninfo->fxn  = (actionFunction) transformCheckOverlap;
  } else if ( actionArgs.CommandIs("contacts")       ) {
    actioninfo->type = TRANSFORM_CONTACTS;
    actioninfo->fxn  = (actionFunction) transformContacts;
  } else if ( actionArgs.CommandIs("correlation")    ) {
    actioninfo->type = TRANSFORM_CORRELATION;
    actioninfo->fxn  = (actionFunction) transformCorr;
  } else if ( actionArgs.CommandIs("clusterdihedral")) {
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
  } else if ( actionArgs.CommandIs("matrix")         ) {
    actioninfo->type = TRANSFORM_MATRIX;
    actioninfo->fxn  = (actionFunction) transformMatrix;
  } else if ( actionArgs.CommandIs("principal")      ) {
    actioninfo->type = TRANSFORM_PRINCIPAL;
    actioninfo->fxn  = (actionFunction) transformPrincipal;
  } else if ( actionArgs.CommandIs("projection")     ) {
    actioninfo->type = TRANSFORM_PROJECTION;
    actioninfo->fxn  = (actionFunction) transformProjection;
  } else if ( actionArgs.CommandIs("randomizeions")  ) {
    actioninfo->type = TRANSFORM_RANDOMIZEIONS;
    actioninfo->fxn  = (actionFunction) transformRandomizeIons;
  } else if ( actionArgs.CommandIs("runningaverage") ) {
    actioninfo->type = TRANSFORM_RUNNINGAVERAGE;
    actioninfo->fxn  = (actionFunction) transformRunningAverage;
  } else if ( actionArgs.CommandIs("scale")          ) {
    actioninfo->type = TRANSFORM_SCALE;
    actioninfo->fxn  = (actionFunction) transformScale;
  } else if ( actionArgs.CommandIs("unwrap")         ) {
    actioninfo->type = TRANSFORM_UNWRAP;
    actioninfo->fxn  = (actionFunction) transformUnwrap;
  } else if ( actionArgs.CommandIs("vector")         ) {
    actioninfo->type = TRANSFORM_VECTOR;
    actioninfo->fxn  = (actionFunction) transformVector;
  } else if ( actionArgs.CommandIs("watershell")     ) {
    actioninfo->type = TRANSFORM_WATERSHELL;
    actioninfo->fxn  = (actionFunction) transformWatershell;
  } else {
    mprinterr("Error: PtrajAction: Unrecognized Ptraj command: %s\n",actionArgs.Command());
    return 1;
  }

  // Convert the ArgList to the argStackType used by functions in ptraj_actions.c
  argumentStack = (argStackType*) malloc( sizeof(argStackType) );
  //mprintf("DEBUG:\targumentStack address (init): %x\n",argumentStack);
  argumentStack->arglist = NULL;
  argumentStack->marked = NULL;
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
  actioninfo->carg1 = (void *) &argumentStack;
  //mprintf("DEBUG:\taction->carg1 address (init): %x\n",actioninfo->carg1);
  // NOTE: Should ptraj be freeing up the args?
  // DEBUG
  //argStackType **argumentStackPointer = (argStackType **) actioninfo->carg1;
  //mprintf("DEBUG:\targumentStack address (init 2): %x\n",*argumentStackPointer);
  if (debug>0) printArgumentStack(&argumentStack);

  // Initialize state memory
  actioninfo->state = (ptrajState*) malloc( sizeof(ptrajState) );
  INITIALIZE_ptrajState( actioninfo->state );

  // Dont call setup here since there is no state information yet
  mprintf("    PTRAJ ACTION: [%s]\n",actionArgs.Command());
  
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
  ptrajState *state = actioninfo->state;
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
  // IPRES - in Cpptraj the atom #s in IPRES (resnums) is shifted by -1
  // to be consistent with the rest of cpptraj; however, ptraj expects
  // standard ipres. Create a copy.
  state->ipres = (int*) malloc( (currentParm->nres+1) * sizeof(int));
  for (int res = 0; res <= currentParm->nres; res++)
    state->ipres[res] = currentParm->resnums[res]+1;
  //state->ipres[currentParm->nres] = currentParm->natom;
  state->ipres_mask = currentParm->resnums;
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
  // Assume if no solvent mask, no solvent present
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
    state->solventAtoms = state->solventMoleculeStop[currentParm->solventMolecules-1] - 
                          state->solventMoleculeStart[0];
  } else {
    memset(state->solventMask,0,currentParm->natom);
    memset(state->solventMoleculeStart,0,currentParm->solventMolecules);
    memset(state->solventMoleculeStop,0,currentParm->solventMolecules);
    state->solventAtoms = 0;
  }
  state->atomName = currentParm->names;
  state->residueName = currentParm->resnames;
  state->maxFrames = DSL->MaxFrames(); // NOTE: Is this correct?
  state->temp0 = 0.0;

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
  x_coord = (double*) malloc( currentParm->natom * sizeof(double));
  y_coord = (double*) malloc( currentParm->natom * sizeof(double));
  z_coord = (double*) malloc( currentParm->natom * sizeof(double));
        
  return 0;  
}

// PtrajAction::action()
int PtrajAction::action() {
  // NOTE: Not checking for Called Setup here since the noSetup
  // flag should catch it.
  // Convert coordinates array X0Y0Z0X1... into separate X Y Z vectors
  int i3 = 0;
  for (int atom = 0; atom < currentFrame->natom; atom++) {
    x_coord[atom] = currentFrame->X[i3++];
    y_coord[atom] = currentFrame->X[i3++];
    z_coord[atom] = currentFrame->X[i3++];
  }
  // Protect state box coords
  ptraj_box[0] = currentFrame->box[0];
  ptraj_box[1] = currentFrame->box[1];
  ptraj_box[2] = currentFrame->box[2];
  ptraj_box[3] = currentFrame->box[3];
  ptraj_box[4] = currentFrame->box[4];
  ptraj_box[5] = currentFrame->box[5];
  
  actioninfo->fxn(actioninfo, x_coord, y_coord, z_coord, ptraj_box, PTRAJ_ACTION);

  return 0;
} 

// PtrajAction::print()
void PtrajAction::print() {
  // NOTE: Is it ok to just call this will NULLs?
  if (CalledSetup)
    actioninfo->fxn(actioninfo, x_coord, y_coord, z_coord, ptraj_box, PTRAJ_PRINT);
}

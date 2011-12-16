#include <cstdlib>
#include <cstring>
#include "ptraj_convert.h"

// CreateArgumentStack()
argStackType *CreateArgumentStack(ArgList &argIn, int debug) {
  argStackType *argumentStack;

  // Convert the ArgList to the argStackType used by functions in ptraj_actions.c
  argumentStack = (argStackType*) malloc( sizeof(argStackType) );
  //mprintf("DEBUG:\targumentStack address (init): %x\n",argumentStack);
  argumentStack->arglist = NULL;
  argumentStack->marked = NULL;
  int nargs = 0;
  char *currentArg = argIn.getNextString();
  while (currentArg != NULL) {
    argumentStack->arglist = (char**) realloc(argumentStack->arglist, (nargs+1) * sizeof(char*));
    argumentStack->arglist[nargs] = (char*) malloc( (strlen(currentArg)+1) * sizeof(char) );
    strcpy(argumentStack->arglist[nargs], currentArg);
    currentArg = argIn.getNextString();
    nargs++;
  } 
  argumentStack->nargs = nargs;
  argumentStack->marked = (char*) malloc(nargs * sizeof(char));
  memset(argumentStack->marked, 'F', nargs);
  if (debug>0) printArgumentStack(&argumentStack);

  return argumentStack;
}

// FreeArgumentStack()
void FreeArgumentStack(argStackType *argumentStack) {
  if (argumentStack!=NULL) {
    // free individual args
    if (argumentStack->arglist!=NULL) {
      for (int arg=0; arg < argumentStack->nargs; arg++)
        if (argumentStack->arglist[arg]!=NULL) free(argumentStack->arglist[arg]);
      free(argumentStack->arglist);
    }
    if (argumentStack->marked!=NULL) free(argumentStack->marked);
    free(argumentStack);
  }
}

// CreateState()
ptrajState *CreateState(AmberParm *currentParm, int maxFrames) {
  ptrajState *state;

  state = (ptrajState*) malloc( sizeof(ptrajState) );
  INITIALIZE_ptrajState( state );

  state->box[0] = currentParm->Box[0];
  state->box[1] = currentParm->Box[1];
  state->box[2] = currentParm->Box[2];
  state->box[3] = currentParm->Box[3];
  state->box[4] = currentParm->Box[4];
  state->box[5] = currentParm->Box[5];
  state->masses = currentParm->mass;
  state->charges = currentParm->Charges_ptr();
  state->atoms = currentParm->natom;
  state->residues = currentParm->Nres();
  // IPRES - in Cpptraj the atom #s in IPRES (resnums) is shifted by -1
  // to be consistent with the rest of cpptraj; however, ptraj actions 
  // expect standard ipres (atom #s start at 1). Create a copy with atom
  // numbers shifted +1.
  int *ResNums = currentParm->ResAtomNums_ptr();
  state->ipres = (int*) malloc( (state->residues+1) * sizeof(int));
  for (int res = 0; res <= state->residues; res++)
    state->ipres[res] = ResNums[res]+1;
  // The Cpptraj version of the mask parser expects Cpptraj-style ipres,
  // so need that as well
  state->ipres_mask = ResNums;
  state->IFBOX = AmberIfbox(currentParm->Box[4]);
  //state->boxfixed
  state->molecules = currentParm->Nmol();
  state->moleculeInfo = currentParm->AtomsPerMol_ptr();
  // Solvent info
  state->solventMask = (int*) malloc( state->atoms * sizeof(int));
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
  state->atomName = currentParm->AtomNames_ptr();
  state->residueName = currentParm->ResidueNames_ptr();
  state->maxFrames = maxFrames; 
  state->temp0 = 0.0;

  return state;
}

// FreeState()
void FreeState(ptrajState *state) {
  if (state!=NULL) {
    if (state->ipres!=NULL) free(state->ipres);
    if (state->solventMask!=NULL) free(state->solventMask);
    if (state->solventMoleculeStart!=NULL) free(state->solventMoleculeStart);
    if (state->solventMoleculeStop!=NULL) free(state->solventMoleculeStop);
    free(state);
  }
}


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
ptrajState *CreateState(Topology *currentParm, int maxFrames) {
  ptrajState *state;

  state = (ptrajState*) malloc( sizeof(ptrajState) );
  INITIALIZE_ptrajState( state );

  // General info
  state->maxFrames = maxFrames; 
  state->temp0 = 0.0;
  currentParm->BoxCoords( state->box );
  if (currentParm->BoxType() == Box::NOBOX)
    state->IFBOX = 0;
  else if (currentParm->BoxIsTruncOct())
    state->IFBOX = 2;
  else
    state->IFBOX = 1;
  //state->boxfixed
  // Atom Info
  state->atoms = currentParm->Natom();
  state->masses = (double*) malloc(state->atoms * sizeof(double)); // ack! malloc!
  state->charges = (double*) malloc(state->atoms * sizeof(double));
  state->atomName = (NAME*) malloc(state->atoms * sizeof(NAME));
  for (int atom = 0; atom < state->atoms; atom++) {
    state->masses[atom] = (*currentParm)[atom].Mass();
    state->charges[atom] = (*currentParm)[atom].Charge();
    (*currentParm)[atom].Name().ToBuffer( state->atomName[atom] );
  }
  // Residue info
  state->residues = currentParm->Nres();
  // IPRES - in Cpptraj the atom #s in IPRES (resnums) is shifted by -1
  // to be consistent with the rest of cpptraj; however, ptraj actions 
  // expect standard ipres (atom #s start at 1). Create a copy with atom
  // numbers shifted +1.
  state->ipres = (int*) malloc( (state->residues+1) * sizeof(int));
  // The Cpptraj version of the mask parser expects Cpptraj-style ipres,
  // so need that as well
  // TODO: Not really needed now, cpptraj has new mask parser
  state->ipres_mask = (int*) malloc( (state->residues+1) * sizeof(int));
  state->residueName = (NAME*) malloc ( (state->residues+1) * sizeof(NAME));
  for (int res = 0; res < state->residues; res++) {
    state->ipres_mask[res] = currentParm->Res(res).FirstAtom();
    state->ipres[res] = state->ipres_mask[res] + 1;
    currentParm->Res(res).Name().ToBuffer( state->residueName[res] );
  }
  state->ipres_mask[state->residues] = state->atoms;
  state->ipres[state->residues] = state->atoms;
  // Molecule info
  state->molecules = currentParm->Nmol();
  state->moleculeInfo = (int*) malloc( state->molecules * sizeof(int));
  int* apm = state->moleculeInfo;
  for (Topology::mol_iterator mol = currentParm->MolStart();
                              mol != currentParm->MolEnd(); mol++)
  {
    *apm = (*mol).NumAtoms();
    ++apm;
  }
  // Solvent info
  state->solventMask = (int*) malloc( state->atoms * sizeof(int));
  for (int atom = 0; atom < state->atoms; atom++)
    state->solventMask[atom] = 0;
  state->solventMolecules = currentParm->Nsolvent();
  state->solventAtoms = 0;
  if (state->solventMolecules > 0) {
    state->solventMoleculeStart = (int*) malloc( state->solventMolecules * sizeof(int));
    state->solventMoleculeStop = (int*) malloc( state->solventMolecules * sizeof(int));
    int smol = 0;
    for (Topology::mol_iterator mol = currentParm->SolventStart();
                                mol != currentParm->SolventEnd(); mol++)
    {
      state->solventMoleculeStart[smol] = (*mol).BeginAtom();
      state->solventMoleculeStop[smol] = (*mol).EndAtom();
      for (int i = state->solventMoleculeStart[smol]; i < state->solventMoleculeStop[smol]; i++) {
        ++state->solventAtoms;
        state->solventMask[i] = 1;
      }
      ++smol;
    }
  }


  return state;
}

// FreeState()
void FreeState(ptrajState *state) {
  if (state!=NULL) {
    if (state->masses!=NULL) free(state->masses);
    if (state->charges!=NULL) free(state->charges);
    if (state->atomName!=NULL) free(state->atomName);
    if (state->ipres!=NULL) free(state->ipres);
    if (state->ipres_mask!=NULL) free(state->ipres_mask);
    if (state->residueName!=NULL) free(state->residueName);
    if (state->moleculeInfo!=NULL) free(state->moleculeInfo);
    if (state->solventMask!=NULL) free(state->solventMask);
    if (state->solventMoleculeStart!=NULL) free(state->solventMoleculeStart);
    if (state->solventMoleculeStop!=NULL) free(state->solventMoleculeStop);
    free(state);
  }
}


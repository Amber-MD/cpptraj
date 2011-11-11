#include <stdio.h>
#include <string.h>
#define PTRAJ_STATE_MODULE
#include "ptraj_state.h"
#undef PTRAJ_STATE_MODULE
#include "ptraj_common.h"
#include "PtrajMask.h"

// SetPrnlev()
/// Set the global ptraj debug level defined in ptraj_arg
// SetPrnlev()
void SetPrnlev(int prnlevIn) {
  prnlev = prnlevIn;
  if (prnlev>0) printf("Info: ptraj prnlev set to %i\n",prnlev);
}

// ========== PTRAJ STATE functions ============================================
ptrajState *ptrajCopyState(ptrajState **stateinp) {
  ptrajState *state, *statein;
  int i;

  /*
   *  Make a copy of the current state
   */
  statein = *stateinp;

  state = (ptrajState *) safe_malloc( sizeof(ptrajState) );
  INITIALIZE_ptrajState(state);

  state->atoms = statein->atoms;
  state->atomName = (Name *) safe_malloc( sizeof(Name) * state->atoms);
  state->masses = (double *) safe_malloc(sizeof(double) * state->atoms);
  state->charges = (double *) safe_malloc(sizeof(double) * state->atoms);
  for (i=0; i < state->atoms; i++) {
    strcpy(state->atomName[i], statein->atomName[i]);
    state->masses[i] = statein->masses[i];
    state->charges[i] = statein->charges[i];
  }

  state->residues = statein->residues;
  state->residueName = (Name *) safe_malloc (sizeof(Name) * state->residues);
  for (i=0; i < state->residues; i++) {
    strcpy(state->residueName[i], statein->residueName[i]);
  }
  state->ipres = (int *) safe_malloc(sizeof(int) * (state->residues+1));
  for (i=0; i <= state->residues; i++) {
    state->ipres[i] = statein->ipres[i];
  }
  state->solventMolecules = statein->solventMolecules;
  state->solventAtoms     = statein->solventAtoms;

  if (statein->solventMolecules) {

    state->solventMask = (int *) safe_malloc(sizeof(int) * statein->atoms);
    for (i=0; i < state->atoms; i++)
      state->solventMask[i] = statein->solventMask[i];

    state->solventMoleculeStart = (int *)
      safe_malloc(sizeof(int) * statein->solventMolecules);
    state->solventMoleculeStop  = (int *)
      safe_malloc(sizeof(int) * statein->solventMolecules);
    for (i=0; i < state->solventMolecules; i++) {
      state->solventMoleculeStart[i] = statein->solventMoleculeStart[i];
      state->solventMoleculeStop[i]  = statein->solventMoleculeStop[i];
    }
  } else {
    state->solventMoleculeStart = NULL;
    state->solventMoleculeStop  = NULL;
    state->solventMask = NULL;
  }

  state->IFBOX = statein->IFBOX;
  for (i=0; i < 6; i++)
    state->box[i] = statein->box[i];

  if ( statein->molecules > 0 ) {

    state->molecules = statein->molecules;
    state->moleculeInfo = (int *)
      safe_malloc( sizeof(int) * state->molecules );
    for (i=0; i < state->molecules; i++) {
      state->moleculeInfo[i] = statein->moleculeInfo[i];
    }
  }

  return(state);
}

// atomToResidue()
int atomToResidue(int atom, int nres, int *resnums) {
  int i;

  if (resnums == NULL) return -1;

  for (i = 0; i < nres; i++)
    if ( atom>=resnums[i] && atom<resnums[i+1] )
      return (i+1);

  return -1;
}


// ========== MASK functions ===================================================

// isActiveDetailed()
int isActiveDetailed(int atom, int residue, int *mask, int atoms, int residues,
                            Name *atomName, Name *residueName, int *ipres)
{
  int i;

  if (residue >= residues || residue < 0) {
    printf("WARNING: residue out of range in isActiveDetailed, res %i (total %i)\n",
           residue, residues);
    return 0;
  }
  for (i = ipres[residue]-1; i < ipres[residue+1]-1; i++)
    if ( mask[i] && strcmp(atomName[i], atomName[atom]) == 0 )
      return 1;

  return 0;

}

// printAtomMaskDetailed()
void printAtomMaskDetailed(FILE *file, int *mask, int atoms, int residues,
                                  Name *atomName, Name *residueName, int *ipres)
{
  //fprintf(stdout,"Warning: printAtomMaskDetailed not implemented for Cpptraj.\n");

  int i, j, curres;
  char tmpatom[20];
  int *resactive, *ressimilar;
  int printed, numactive, numresactive, numressimilar;
  int incurres, innextres;

  printed = 0;
  numactive = 0;
  numresactive = 0;
  numressimilar = 0;

  /*
   *  This routine is kind of junky since in general I want to avoid printing
   *  as much detail as possible.  Therefore, we check to see is certain ranges 
   *  residues have all atoms active, etc. to avoid printing each atom in a residue.
   *  This makes it ugly and obtuse.
   */


  if (mask == NULL) {
    fprintf(file, "[No atoms are selected]");
    return;
  }
 j=0;
  for (i=0; i < atoms; i++)
    if (mask[i]) j++;

  if (j == 0) {
    fprintf(file, "[No atoms are selected]");
    return;
  }

     /*
      *  check if all atoms are active and if so print an asterisk
      */

  j = 0;
  for (i=0; i < atoms; i++) {
    j += mask[i];
  }
  if ( j == atoms ) {
    fprintf(file, "  * (All atoms are selected)");
    return;
  }
  numactive = j;

     /*
      *  determine which residues have all the atoms in that residue active
      */

  resactive = (int *) safe_malloc(sizeof(int) * residues);
  for (i=0; i < residues; i++) {
    resactive[i] = 0.0;
  }
 curres = 0;
  j = 0;
  for (i=0; i < atoms; i++) {
    if (i == ipres[curres+1]-1) {
      if (j == ipres[curres+1] - ipres[curres]) {
        resactive[curres] = 1.0;
        numresactive++;
      }
      j = 0;
      curres++;
    }
    if (mask[i])
      j++;
  }
  //fprintf(stdout,"DEBUG:\tprintAtomMaskDetailed: curres = %i\n",curres);
  // NOTE: DRR - unneccessary, curres will end up == residues
  //if (j == ipres[curres+1] - ipres[curres]) {
  //  resactive[curres] = 1.0;
  //}

     /*
      *  determine the range over which the residues are fully active
      */

  for (curres = residues-2; curres >= 0; curres--) {
    if (resactive[curres]) {
      resactive[curres] += resactive[curres+1];
      numresactive--;
    }
  }


 /*
   *  determine ranges over which residues have the same atoms active
   *  as the next residue
   */
  ressimilar = (int *) safe_malloc(sizeof(int) * residues);
  for (i=0; i < residues; i++) {
    ressimilar[i] = 0.0;
  }

  for (curres = residues-2; curres >=0; curres--) {

    incurres = 0;
    innextres = 0;
    for (i = ipres[curres]-1; i < ipres[curres+2]-1; i++) { /* check current and next residue */
      if ( mask[i] ) {
        incurres++;
      }

      if (isActiveDetailed(i, i<ipres[curres+1]-1 ? curres+1 : curres, mask,
                           atoms, residues, atomName, residueName, ipres))
        innextres++;

      if (incurres != innextres) /* select counterparts in next residues too! */
        break;
    }
    if (incurres && innextres == incurres) {
      ressimilar[curres] = ressimilar[curres+1] + 1;
    } else {
      numressimilar++;
    }

  }
   /*
      *  do the actual printing
      */

  j = 0;
  for (curres = 0; curres < residues; curres++) {

    if (resactive[curres] ) {

      /*
       *  If all of the atoms are active in this residue print either the
       *  residue number or range as appropriate
       */

      if (resactive[curres] > 2) {
        if (j!=0 && j%10 != 0) fprintf(file, ",");
        fprintf(file, ":%i-%i", curres+1, curres+resactive[curres]);
        curres += resactive[curres]-1;
      } else {
        if (j!=0 && j%10 != 0) fprintf(file, ",");
        fprintf(file, ":%i", curres+1);
      }
      j++;
      if (j != 0 && j % 10 == 0) {
        fprintf(file, "\n    ");
        j = 0;
      }
    } else if (ressimilar[curres]) {

      /*
       *  If there is a set of residues with a similar atom selection...
       */
      printed = 0;
      if (ressimilar[curres] >= 1) {
        if (j!=0 && j%10 != 0) fprintf(file, ",");
        fprintf(file, ":%i-%i", curres+1, curres+ressimilar[curres]+1);
        curres += ressimilar[curres];
      } else {
        if (j!=0 && j%10 != 0) fprintf(file, ",");
        fprintf(file, ":%i", curres+1);
      }

      for (i = ipres[curres]-1; i < ipres[curres+1]-1; i++) {
        if ( mask[i] ) {
          if (printed)
            fprintf(file, ",");
          else {
            fprintf(file, "@");
            printed = 1;
          }
          strcpy(tmpatom, atomName[i]);
          fprintf(file, "%s", strtok(tmpatom, " "));
        }
      }
      j++;
      if (j != 0 && j % 10 == 0) {
        fprintf(file, "\n    ");
        j = 0;
      }


    } else {

      /*
       *  Print individual atoms
       */
      if (numactive > 10 && numressimilar > 10 && numresactive > 10) {
        fprintf(file, "\n    ");
        numactive = 0;
      }
      for (i = ipres[curres]-1; i < ipres[curres+1]-1; i++) {
        if ( mask[i] ) {
          if (j!=0 && j%10 != 0) fprintf(file, ",");

          strcpy(tmpatom, atomName[i]);
          fprintf(file, ":%i@%s", curres+1, strtok(tmpatom, " "));
          j++;
        }
        if (j != 0 && j % 10 == 0) {
          fprintf(file, "\n    ");
          j = 0;
        }
      }
    }
  }
  /*
    if (j!=0 && j%10 != 0) fprintf(file, "\n");
  */
  safe_free(resactive);
  safe_free(ressimilar);
}

// printAtomMask()
void printAtomMask(FILE *file, int *mask, ptrajState *state) {
  printAtomMaskDetailed(file, mask, state->atoms, state->residues, state->atomName,
                        state->residueName, state->ipres);
}

// printAtomCompact2()
void printAtomCompact2(FILE *fpout, int atom, ptrajState *state) {
  int curres;
  char buffer[50], *bufferp;

  curres = atomToResidue(atom+1, state->residues, state->ipres)-1;

  sprintf(buffer, ":%i", curres+1);
  bufferp = buffer+strlen(buffer);
  sprintf(bufferp, "@%s", state->atomName[atom]);

  fprintf(fpout, "%s", buffer);
}


// processAtomMask()
int *processAtomMask( char *maskString, ptrajState *state ) {
  char *charmask;
  int *intmask;
  int i, actualAtoms;

  intmask = (int *) safe_malloc( state->atoms * sizeof(int) );
  
  if ( strcmp(maskString, "") == 0) {
    maskString = strcpy(maskString, "*");
  }

 charmask = parseMaskC(maskString, state->atoms, state->residues,
                        state->atomName, state->residueName, state->ipres_mask,
                        NULL, NULL, prnlev);
  //  The new mask parsing routine returns a character array
  //  rather than integer; this makes more sense, however in the meantime
  //  we need to convert between the two.
  actualAtoms = 0;
  for (i = 0; i < state->atoms; i++) {
    if (charmask[i] == 'T') { 
      intmask[i] = 1;
      actualAtoms++;
    } else
      intmask[i] = 0;
  }
  
  if (actualAtoms > 0) {
    fprintf(stdout, "Mask [%s] represents %i atoms\n",
            maskString, actualAtoms);
  } else  {
    fprintf(stdout, "Mask [%s] represents %i atoms ",
            maskString, actualAtoms);
    fprintf(stdout, "!!!NO ATOMS DETECTED!!!\n");
    safe_free(intmask);
    intmask = NULL;
  }
  // DEBUG
  //fprintf(stdout,"Mask:\n");
  //for (i = 0; i < state->atoms; i++) 
  //  fprintf(stdout,"\t%8i[%s]: %i\n",i,state->atomName[i],intmask[i]);

  if (prnlev > 2) {
    fprintf(stdout, "Parsed mask string matches:\n");
    printAtomMaskDetailed(stdout, intmask, state->atoms, state->residues, state->atomName, 
                          state->residueName, state->ipres);
  }

  safe_free(charmask);
  return(intmask);
}


#include <cstdlib>
#include <cstring>
#include "AtomMask.h"
#include "CpptrajStdio.h"
#include "ptrajmask.h"

// CONSTRUCTOR
AtomMask::AtomMask() {
  invertMask=false;
  maskString=NULL;
  Selected=NULL;
  Nselected=0;
  CharMask=NULL;
}

// DESTRUCTOR
AtomMask::~AtomMask() {
  if (maskString!=NULL) free(maskString);
  if (Selected!=NULL) free(Selected);
  if (CharMask!=NULL) free(CharMask);
}

/*
 * AtomMask::AddAtom()
 * Add atom to Selected array in this mask.
 */
void AtomMask::AddAtom(int atom) {
  Selected = (int*) realloc(Selected, (Nselected+1) * sizeof(int));
  Selected[Nselected] = atom;
  Nselected++;
}

/*
 * AtomMask::PrintMaskAtoms()
 * Print all atoms in mask to line.
 */
void AtomMask::PrintMaskAtoms() {
  if (this->None()) 
    mprintf("No atoms selected.");
  else {
    for (int atom=0; atom<Nselected; atom++)
      mprintf(" %i",Selected[atom]);
  }
}

/*
 * AtomMask::Copy()
 * Return a copy of this atom mask.
 */
AtomMask *AtomMask::Copy() {
  int mask;
  AtomMask *newMask;
  
  newMask = new AtomMask();
  if (this->Selected!=NULL) {
    newMask->Selected = (int*) malloc(this->Nselected * sizeof(int));
    for (mask=0; mask < this->Nselected; mask++) 
      newMask->Selected[mask] = this->Selected[mask];
    newMask->Nselected = this->Nselected;
  }
  newMask->invertMask = this->invertMask;
  if (this->maskString!=NULL) {
    newMask->maskString = (char*) malloc( (strlen(this->maskString)+1) * sizeof(char));
    strcpy(newMask->maskString, this->maskString);
  }
  if (this->CharMask!=NULL) {
    // Nselected is set to natom, the size of char mask array
    newMask->CharMask = (char*) malloc(this->Nselected * sizeof(char));
    strcpy(newMask->CharMask, this->CharMask);
    newMask->Nselected = this->Nselected;
  }

  return newMask;
}

/*
 * AtomMask::SetMaskString()
 * Set maskString, replacing any existing maskString 
 * If maskStringIn is NULL set to * (all atoms)
 */
void AtomMask::SetMaskString(char *maskStringIn) {
  if (maskString!=NULL) free(maskString);
  if (maskStringIn!=NULL) {
    maskString = (char*) malloc( (strlen(maskStringIn)+1) * sizeof(char));
    strcpy(maskString, maskStringIn);
  } else {
    maskString = (char*) malloc( 2 * sizeof(char));
    strcpy(maskString, "*");
  }
}

/*
 * AtomMask::None()
 * Return true if no atoms are selected.
 */
bool AtomMask::None() {
  if (Nselected==0) return true;
  return false;
}

/*
 * AtomMask::SetupMask()
 * Set up an atom mask given a parm file. The basic atom mask is allocated
 * using PTRAJs old mask parser (accessed from the AmberParm) which returns
 * a char array, 1 for each atom, where selected atoms are denoted by T.
 * Base on this create an array of selected atom numbers. If inverMask is true 
 * create an array of atoms that are not selected.
 */
int AtomMask::SetupMask(AmberParm *Pin, int debug) {
  int atom;
  char *mask;
  char maskChar;

  if (Pin==NULL) {
    mprintf("    Error: AtomMask::SetupMask: (%s) Topology is NULL.\n", maskString);
    return 1;
  }
  maskChar='T';
  if (invertMask) maskChar='F';

  // Allocate atom mask - free mask if already allocated
  // NOTE: Args 7-8 are for distance criteria selection. Arg 7 is coords,
  //       arg 8 should be f for float (NOT IMPLEMENTED) or d for double.
  //       Last arg is debug level.
  mask = parseMaskString(maskString, Pin->natom, Pin->nres, Pin->names, Pin->resnames,
                         Pin->resnums, NULL, Pin->types, debug);
  if (mask==NULL) {
    mprintf("    Error: Could not set up mask %s for topology %s\n",
            maskString, Pin->parmName);
    return 1;
  }

  // Set up an integer list of the selected atoms. 
  // NOTE: For large selections this will use 4x the memory of the char atom
  //       mask. Could check to see which will be bigger.
  Nselected=0;
  if (Selected!=NULL) free(Selected);
  Selected = (int*) malloc( Pin->natom * sizeof(int) );
  for (atom=0; atom<Pin->natom; atom++) {
    if (mask[atom]==maskChar) {
      Selected[Nselected]=atom;
      Nselected++;
    }
  }
  // Resize array for number of selected atoms
  Selected = (int*) realloc(Selected, Nselected * sizeof(int) );

  if (debug>0) {
    if (invertMask)
      mprintf("          Inverse of Mask %s corresponds to %i atoms.\n",
              maskString,Pin->natom - Nselected);
    else
      mprintf("          Mask %s corresponds to %i atoms.\n",maskString,Nselected);
  }

  // Free the character mask, no longer needed
  free(mask);
  
  return 0;
}

/*
 * AtomMask::SetupCharMask()
 * For cases where we need to know both atoms in and out of mask
 * just use the old school char array. In this case Nselected will
 * be the total size of the mask, not just the # of selected atoms.
 */
int AtomMask::SetupCharMask(AmberParm *Pin, int debug) {
  if (Pin==NULL) {
    mprintf("    Error: AtomMask::SetupCharMask: (%s) Topology is NULL.\n", maskString);
    return 1;
  }

  // Allocate atom mask - free mask if already allocated
  Nselected = 0;
  if (CharMask!=NULL) free(CharMask);
  CharMask = parseMaskString(maskString, Pin->natom, Pin->nres, Pin->names, Pin->resnames,
                             Pin->resnums, NULL, Pin->types, debug);
  if (CharMask==NULL) {
    mprintf("    Error: Could not set up mask %s for topology %s\n",
            maskString, Pin->parmName);
    return 1;
  }
  Nselected = Pin->natom;
  return 0;
}

/*
 * AtomMask::SetupCharMask()
 * Set up old school char array with coordinates.
 */
int AtomMask::SetupCharMask(AmberParm *Pin, double *Xin, int debug) {
  if (Pin==NULL) {
    mprintf("    Error: AtomMask::SetupCharMask: (%s) Topology is NULL.\n", maskString);
    return 1;
  }

  // Allocate atom mask - free mask if already allocated
  Nselected = 0;
  if (CharMask!=NULL) free(CharMask);
  CharMask = parseMaskString(maskString, Pin->natom, Pin->nres, Pin->names, Pin->resnames,
                             Pin->resnums, Xin, Pin->types, debug);
  if (CharMask==NULL) {
    mprintf("    Error: Could not set up mask %s for topology %s\n",
            maskString, Pin->parmName);
    return 1;
  }
  Nselected = Pin->natom;
  return 0;
}

/*
 * AtomMask::AtomInCharMask()
 * If CharMask has been set up check if atom has been selected.
 */
bool AtomMask::AtomInCharMask(int atom) {
  if (CharMask==NULL) return false;
  if (atom < 0) return false;
  if (atom >= Nselected) return false;
  if (CharMask[atom]=='T') return true;
  return false;
}

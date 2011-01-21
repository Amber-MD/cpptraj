#include <cstdlib>
#include <cstring>
#include "AtomMask.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
AtomMask::AtomMask() {
  invertMask=false;
  maskString=NULL;
  Selected=NULL;
  Nselected=0;
  //P=NULL;
  //N=0;
  debug=0;
}

// DEBUG CONSTRUCTOR
/*
AtomMask::AtomMask(int debugIn) {
  invertMask=false;
  maskString=NULL;
  Selected=NULL;
  Nselected=0;
  P=NULL;
  N=0;
  debug=debugIn;
}*/

// DESTRUCTOR
AtomMask::~AtomMask() {
  if (maskString!=NULL) free(maskString);
  if (Selected!=NULL) free(Selected);
}

/*
 * AtomMask::AddAtom()
 * Add atom to this mask.
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
 */
AtomMask *AtomMask::Copy() {
  int mask;
  AtomMask *newMask;
  
  newMask = new AtomMask();
  newMask->Selected = (int*) malloc(this->Nselected * sizeof(int));
  for (mask=0; mask < this->Nselected; mask++) 
    newMask->Selected[mask] = this->Selected[mask];
  newMask->Nselected = this->Nselected;
  newMask->invertMask = this->invertMask;
  if (this->maskString!=NULL) {
    newMask->maskString = (char*) malloc( (strlen(this->maskString)+1) * sizeof(char));
    strcpy(newMask->maskString, this->maskString);
  }
  //newMask->P = this->P;
  //newMask->N = this->N;
  newMask->debug = this->debug;

  return newMask;
}

/*
 * AtomMask::SetMaskString()
 * Set maskString, replacing any existing maskString and Selection. 
 * If maskStringIn is NULL set to * (all atoms)
 */
void AtomMask::SetMaskString(char *maskStringIn) {
  if (maskString!=NULL) free(maskString);
  if (Selected!=NULL) free(Selected);
  Selected=NULL;
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
 * NOTE: Do we really need to store the Parm?
 */
int AtomMask::SetupMask(AmberParm *Pin, int debugIn) {
  int atom;
  //size_t SelectedSize, maskSize;
  char *mask;
  char maskChar;

  if (Pin==NULL) {
    mprintf("    Error: AtomMask::SetupMask: (%s) Topology is NULL.\n", maskString);
    return 1;
  }
  debug=debugIn;
  maskChar='T';
  if (invertMask) maskChar='F';

  // Allocate atom mask - free mask if already allocated
  //P = Pin;
  mask = Pin->mask(maskString);
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

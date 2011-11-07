#include "AtomMask.h"
#include "CpptrajStdio.h"
#include "PtrajMask.h"
#include <cstring>
#include <cstdlib> // needed for malloc / free, ptraj_actions.c

// CONSTRUCTOR
AtomMask::AtomMask() {
  invertMask=false;
  maskString=NULL;
  Selected=NULL;
  Nselected=0;
  Nchar=0;
  CharMask=NULL;
  Postfix = NULL;
}

// DESTRUCTOR
AtomMask::~AtomMask() {
  if (maskString!=NULL) delete[] maskString;
  if (Selected!=NULL) delete[] Selected;
  //if (CharMask!=NULL) delete[] CharMask;
  if (CharMask!=NULL) free(CharMask);
  if (Postfix!=NULL) delete[] Postfix;
}

/* AtomMask::Reset()
 */
void AtomMask::Reset() {
  if (maskString!=NULL) delete[] maskString;
  if (Selected!=NULL) delete[] Selected;
  //if (CharMask!=NULL) delete[] CharMask;
  if (CharMask!=NULL) free(CharMask);
  invertMask=false;
  maskString=NULL;
  Selected=NULL;
  Nselected=0;
  Nchar=0;
  CharMask=NULL;
}

/* AtomMask::AddAtom()
 * Add atom to Selected array in this mask.
 */
void AtomMask::AddAtom(int atom) {
  // Ensure atom is not already in mask
  for (int maskidx=0; maskidx < Nselected; maskidx++)
    if (Selected[maskidx]==atom) return;
  int newN = Nselected + 1;
  int *newSelected = new int[ newN ];
  if (Selected!=NULL) {
    memcpy(newSelected,Selected,Nselected*sizeof(int));
    delete[] Selected;
  }
  newSelected[Nselected] = atom;
  Selected = newSelected;
  Nselected = newN;
}

/* AtomMask::AddAtoms()
 * Given an array, add the atom numbers in array to the Selected array.
 * NOTE: Unlike AddAtom the incoming list is not checked for duplicates.
 */
void AtomMask::AddAtoms(int *atomList, int N_to_add) {
  if (atomList==NULL || N_to_add < 1) return;
  // Allocate space for new array
  int newN = Nselected + N_to_add;
  int *newSelected = new int[ newN ];
  // Copy over old selected array
  if (Selected!=NULL) {
    memcpy(newSelected, Selected, Nselected * sizeof(int));
    delete[] Selected;
  }
  // Add atoms in atomList
  for (int atom = 0; atom < N_to_add; atom++) 
    newSelected[Nselected++] = atomList[atom];
  Selected = newSelected;
  Nselected = newN;
}

/* AtomMask::AddAtomRange()
 * Add atoms in range from minAtom up to but not including maxAtom to 
 * Selected array.
 * NOTE: Does not check for duplicates.
 */
void AtomMask::AddAtomRange(int minAtom, int maxAtom) {
  if (minAtom >= maxAtom) return;
  int newN = Nselected + (maxAtom - minAtom);
  int *newSelected = new int[ newN ];
  // Copy over old selected array
  if (Selected!=NULL) {
    memcpy(newSelected, Selected, Nselected * sizeof(int));
    delete[] Selected;
  }
  // Add minAtom <= atom < maxAtom to Selected
  for (int atom = minAtom; atom < maxAtom; atom++) 
    newSelected[Nselected++] = atom;
  Selected = newSelected;
  Nselected = newN;
}

/* AtomMask::PrintMaskAtoms()
 * Print all atoms in mask to line.
 */
void AtomMask::PrintMaskAtoms() {
  if (this->None()) 
    mprintf("No atoms selected.");
  else if (Selected!=NULL) {
    for (int atom=0; atom<Nselected; atom++)
      mprintf(" %i",Selected[atom]);
  } else if (CharMask!=NULL) {
    for (int atom=0; atom<Nchar; atom++)
      if (CharMask[atom]=='T') mprintf(" %i",atom);
  } else 
    mprintf("Warning: Mask [%s] has not been set up yet.\n",maskString);
}

/* AtomMask::operator=()
 */
AtomMask &AtomMask::operator=(const AtomMask &rhs) {
  // Check for self assignment
  if ( this == &rhs ) return *this;

  // Deallocate
  //if (CharMask!=NULL) delete[] CharMask;
  if (CharMask!=NULL) free(CharMask);
  if (maskString!=NULL) delete[] maskString;
  if (Selected!=NULL) delete[] Selected;

  // Allocate and copy
  invertMask = rhs.invertMask;
  Nchar = rhs.Nchar;
  Nselected = rhs.Nselected;
  if (rhs.CharMask!=NULL) {
    //CharMask = new char[ Nchar ];
    CharMask = (char*) malloc( Nchar * sizeof(char));
    memcpy(CharMask, rhs.CharMask, Nchar * sizeof(char));
  }
  if (rhs.maskString!=NULL) {
    maskString = new char[ strlen(rhs.maskString) + 1 ];
    strcpy(maskString, rhs.maskString);
  }
  if (rhs.Selected!=NULL) {
    Selected = new int[ Nselected ];
    memcpy(Selected, rhs.Selected, Nselected * sizeof(int));
  }

  // Return *this
  return *this;
}

/* AtomMask::CopyMask()
 * Return a copy of this atom mask.
 */
AtomMask *AtomMask::CopyMask() {
  AtomMask *newMask;
  
  newMask = new AtomMask();

  if (this->Selected!=NULL) {
    newMask->Selected = new int[ this->Nselected ];
    memcpy(newMask->Selected, this->Selected, this->Nselected * sizeof(int));
    newMask->Nselected = this->Nselected;
  }

  newMask->invertMask = this->invertMask;

  if (this->maskString!=NULL) {
    newMask->maskString = new char[ strlen(this->maskString)+1 ];
    strcpy(newMask->maskString, this->maskString);
  }

  if (this->CharMask!=NULL) {
    //newMask->CharMask = new char[ this->Nchar ];
    newMask->CharMask = (char*) malloc(this->Nchar * sizeof(char));
    memcpy(newMask->CharMask, this->CharMask, this->Nchar * sizeof(char));
    newMask->Nchar = this->Nchar;
    newMask->Nselected = Nselected;
  }

  return newMask;
}

/* AtomMask::SetMaskString()
 * Set maskString, replacing any existing maskString 
 * If maskStringIn is NULL set to * (all atoms)
 */
int AtomMask::SetMaskString(char *maskStringIn) {
  char infix[MAXSELE], postfix[MAXSELE];
  int debug = 0; // temporary hack

  if (maskString!=NULL) delete[] maskString;
  if (Postfix!=NULL) delete[] Postfix;
  Postfix = NULL;
  if (maskStringIn!=NULL) {
    maskString = new char[ strlen(maskStringIn)+1 ];
    strcpy(maskString, maskStringIn);
  } else {
    maskString = new char[2];
    maskString[0] = '*';
    maskString[1] = '\0'; 
  }

  // 1) preprocess input expression
  if (tokenize(maskString, infix)!=0) {
    delete[] maskString; 
    maskString=NULL;
    return 1;
  }
  if (debug > 5) mprintf("tokenized: ==%s==\n", infix);

  // 2) construct postfix (RPN) notation 
  if (torpn(infix, postfix)!=0) {
    delete[] maskString;
    maskString=NULL;
    return 1;
  }
  if (debug > 5) mprintf("postfix  : ==%s==\n", postfix);
  Postfix = new char[ strlen(postfix)+1 ];
  strcpy(Postfix,postfix);
  return 0;
}

/* AtomMask::None()
 * Return true if no atoms are selected.
 */
bool AtomMask::None() {
  if (Nselected==0) return true;
  return false;
}

/* AtomMask::SetupMask()
 * Set up an atom mask given a parm file. The basic atom mask is allocated
 * using a version of PTRAJs mask parser (PtrajMask::parseMaskString)
 * which returns a char array of size P->natom, where selected atoms are 
 * denoted by T. Based on this create an array of selected atom numbers. 
 * If invertMask is true create an array of atoms that are not selected
 * instead.
 */
int AtomMask::SetupMask(AmberParm *Pin, double *Xin, int debug) {
  char *mask;
  char maskChar = 'T';

  if (Pin==NULL) {
    mprinterr("    Error: AtomMask::SetupMask: (%s) Topology is NULL.\n", maskString);
    return 1;
  }
  if (Postfix==NULL) {
    mprinterr("    Error: AtomMask::SetupMask: Postfix is NULL.\n");
    return 1;
  }

  if (invertMask) maskChar='F';

  if (debug>1) mprintf("\tAtomMask: Parsing postfix [%s]\n",Postfix);

  // Allocate character atom mask
  mask = parseMaskString(Postfix, Pin->natom, Pin->nres, Pin->names, Pin->resnames,
                         Pin->resnums, Xin, Pin->types, debug);
  if (mask==NULL) {
    mprinterr("    Error: Could not set up mask %s for topology %s\n",
            maskString, Pin->parmName);
    return 1;
  }

  // Wipe out previous char mask if allocated
  Nchar = 0;
  //if (CharMask!=NULL) delete[] CharMask;
  if (CharMask!=NULL) free(CharMask);
  CharMask = NULL;
  // Set up an integer list of the selected atoms. 
  // NOTE: For large selections this will use 4x the memory of the char atom
  //       mask. Could check to see which will be bigger.
  Nselected=0;
  if (Selected!=NULL) delete[] Selected;
  Selected = new int[ Pin->natom ];
  for (int atom=0; atom<Pin->natom; atom++) {
    if (mask[atom]==maskChar) {
      Selected[Nselected]=atom;
      Nselected++;
    }
  }
  // Resize array for number of selected atoms
  int *newSelected = new int[ Nselected ];
  memcpy(newSelected, Selected, Nselected * sizeof(int));
  delete[] Selected;
  Selected = newSelected;

  if (debug>0) {
    if (invertMask)
      mprintf("          Inverse of Mask %s corresponds to %i atoms.\n",
              maskString,Pin->natom - Nselected);
    else
      mprintf("          Mask %s corresponds to %i atoms.\n",maskString,Nselected);
  }

  // Free the character mask, no longer needed.
  // NOTE: Use free() for now while PtrajMask uses free, needed by
  //       ptraj_actions.c
  //delete[] mask;
  free(mask);
  
  return 0;
}

/* AtomMask::SetupCharMask()
 * For cases where we need to know both atoms in and out of mask
 * just use the old school char array. 
 */
/*int AtomMask::SetupCharMask(AmberParm *Pin, int debug) {
  if (Pin==NULL) {
    mprintf("    Error: AtomMask::SetupCharMask: (%s) Topology is NULL.\n", maskString);
    return 1;
  }

  // Allocate atom mask - free mask if already allocated
  Nchar = 0;
  if (CharMask!=NULL) delete[] CharMask;
  CharMask = parseMaskString(Postfix, Pin->natom, Pin->nres, Pin->names, Pin->resnames,
                             Pin->resnums, NULL, Pin->types, debug);
  if (CharMask==NULL) {
    mprintf("    Error: Could not set up mask %s for topology %s\n",
            maskString, Pin->parmName);
    return 1;
  }
  // Determine number of selected atoms
  Nselected=0;
  for (int atom=0; atom<Pin->natom; atom++)
    if (CharMask[atom]=='T') Nselected++;
  Nchar = Pin->natom;
  return 0;
}*/

/* AtomMask::SetupCharMask()
 * For cases where we need to know both atoms in and out of mask use
 * the output of parseMaskString (1 char for each atom, T for selected,
 * F for not.). 
 */
int AtomMask::SetupCharMask(AmberParm *Pin, double *Xin, int debug) {
  if (Pin==NULL) {
    mprinterr("    Error: AtomMask::SetupCharMask: (%s) Topology is NULL.\n", maskString);
    return 1;
  }
  if (Postfix==NULL) {
    mprinterr("    Error: AtomMask::SetupCharMask: Postfix is NULL.\n");
    return 1;
  }

  // Wipe out previous Selected mask if allocated
  if (Selected!=NULL) delete[] Selected;
  Selected = NULL;

  // Allocate atom mask - free mask if already allocated
  Nchar = 0;
  //if (CharMask!=NULL) delete[] CharMask;
  if (CharMask!=NULL) free( CharMask);
  CharMask = parseMaskString(Postfix, Pin->natom, Pin->nres, Pin->names, Pin->resnames,
                             Pin->resnums, Xin, Pin->types, debug);
  if (CharMask==NULL) {
    mprinterr("    Error: Could not set up mask %s for topology %s\n",
            maskString, Pin->parmName);
    return 1;
  }
  // Determine number of selected atoms
  Nselected=0;
  for (int atom=0; atom<Pin->natom; atom++)
    if (CharMask[atom]=='T') Nselected++;
  Nchar = Pin->natom;
  return 0;
}

/* AtomMask::AtomInCharMask()
 * If CharMask has been set up check if atom has been selected.
 */
bool AtomMask::AtomInCharMask(int atom) {
  if (CharMask==NULL) return false;
  if (atom < 0) return false;
  if (atom >= Nchar) return false;
  if (CharMask[atom]=='T') return true;
  return false;
}

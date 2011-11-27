#include <algorithm> // sort, unique
#include "AtomMask.h"
#include "CpptrajStdio.h"
#include "PtrajMask.h"

// CONSTRUCTOR
AtomMask::AtomMask() {
  Nselected=0;
  maskChar = 'T';
}

// DESTRUCTOR
AtomMask::~AtomMask() {
}

// COPY CONSTRUCTOR
AtomMask::AtomMask(const AtomMask &rhs) {
  Nselected=rhs.Nselected;
  maskChar=rhs.maskChar;
  maskString=rhs.maskString;
  Selected=rhs.Selected;
  CharMask=rhs.CharMask;
  Postfix=rhs.Postfix;
}

// AtomMask::operator=()
/// AtomMask assignment
AtomMask &AtomMask::operator=(const AtomMask &rhs) {
  // Check for self assignment
  if ( this == &rhs ) return *this;
  // Deallocate
  // Allocate and copy
  Nselected=rhs.Nselected;
  maskChar=rhs.maskChar;
  maskString=rhs.maskString;
  Selected=rhs.Selected;
  CharMask=rhs.CharMask;
  Postfix=rhs.Postfix;
  // Return *this
  return *this;
}

// AtomMask::ResetMask()
void AtomMask::ResetMask() {
  Nselected = 0;
  maskChar = 'T';
  maskString.clear();
  Selected.clear();
  CharMask.clear();
  Postfix.clear();
}

// AtomMask::InvertMask()
/** Currently for integer masks only. Reverse the selected mask char for 
  * next selection. By default atoms in the Selected array are those marked 
  * by 'T'. After a call to InvertMask the atoms in the Selected array will
  * be those marked by 'F'. Subsqeuent calls flip the character back and 
  * forth.
  */
void AtomMask::InvertMask() {
  if (maskChar == 'T')
    maskChar = 'F';
  else
    maskChar = 'T';
}

// AtomMask::CopyMask()
AtomMask *AtomMask::CopyMask() {
  AtomMask *newMask;
  newMask = new AtomMask();
  newMask->Nselected = Nselected;
  newMask->maskChar = maskChar;
  newMask->maskString = maskString;
  newMask->Selected = Selected;
  newMask->CharMask = CharMask;
  newMask->Postfix = Postfix;
  return newMask;
}

// AtomMask::AddAtom()
/** Attempt to enforce some sorting by looking for the atom in the mask;
  * as soon as an atom # is found larger than atomIn, insert it at the
  * previous spot.
  */
// 32 33 34 48 49 50
void AtomMask::AddAtom(int atomIn) {
  // Ensure atom is not already in mask
  for (std::vector<int>::iterator atom = Selected.begin();  atom != Selected.end(); atom++) {
    if ( *atom == atomIn) return;
    if ( *atom > atomIn) {
      // Insert at the current position, which is the first atom # > atomIn
      Selected.insert(atom, atomIn);
      return;
    }
  }

  // Add atom to mask
  Selected.push_back(atomIn);
}

// AtomMask::AddAtoms()
/** Given an array, add the atom numbers in array to the Selected array.
  * The resulting array is sorted and any duplicates are removed.
  */
void AtomMask::AddAtoms(std::vector<int> &atomsIn) {
  std::vector<int>::iterator atom;
  // Make room for atomsIn in Selected
  //Selected.reserve( Selected.size() + atomsIn.size() );
  // Put every atom in atomsIn in Selected array
  for (atom = atomsIn.begin(); atom != atomsIn.end(); atom++) 
    Selected.push_back( *atom ); 
  // Sort Selected
  sort( Selected.begin(), Selected.end() );
  // Remove duplicates
  atom = unique( Selected.begin(), Selected.end() );
  Selected.resize( atom - Selected.begin() );
}

// AtomMask::AddAtomRange()
/** Add atoms in range from minAtom up to but not including maxAtom to 
  * Selected array. The resulting array is sorted and duplicates are removed.
  */
void AtomMask::AddAtomRange(int minAtom, int maxAtom) {
  if (minAtom >= maxAtom) return;
  for (int atom = minAtom; atom < maxAtom; atom++)
    Selected.push_back( atom );
  // Sort Selected
  sort( Selected.begin(), Selected.end() );
  // Remove duplicates
  std::vector<int>::iterator atomit = unique( Selected.begin(), Selected.end() );
  Selected.resize( atomit - Selected.begin() );
}

// AtomMask::PrintMaskAtoms()
void AtomMask::PrintMaskAtoms() {
  if (this->None()) 
    mprintf("No atoms selected.");
  else if (!Selected.empty()) {
    for (std::vector<int>::iterator atom = Selected.begin();  atom != Selected.end(); atom++)
      mprintf(" %i",*atom);
  } else if (!CharMask.empty()) {
    int atomnum = 0;
    for (std::vector<char>::iterator mc = CharMask.begin();  mc != CharMask.end(); mc++) {
      if (*mc == maskChar) mprintf(" %i",atomnum);
      ++atomnum;
    }
  } else 
    mprintf("Warning: Mask [%s] has not been set up yet.\n",maskString.c_str());
}

// AtomMask::SetMaskString()
/** Take the given mask expression and preprocess it for subsequent use
  * with the mask parser. Convert to infix, then postfix notation.
  */
int AtomMask::SetMaskString(char *maskStringIn) {
  char infix[MAXSELE], postfix[MAXSELE];
  int debug = 0; // temporary hack

  if (maskStringIn!=NULL) 
    maskString.assign( maskStringIn );
  else
    maskString.assign( "*" );

  // 1) preprocess input expression
  if (tokenize((char*)maskString.c_str(), infix)!=0) 
    return 1;
  
  if (debug > 5) mprintf("tokenized: ==%s==\n", infix);

  // 2) construct postfix (RPN) notation 
  if (torpn(infix, postfix)!=0) 
    return 1;
  
  if (debug > 5) mprintf("postfix  : ==%s==\n", postfix);
  Postfix.assign( postfix );
  return 0;
}

// AtomMask::None()
bool AtomMask::None() {
  if (Nselected==0) return true;
  return false;
}

// AtomMask::SetupMask()
/** Set up an atom mask given a parm file. The basic atom mask is allocated
  * using a version of PTRAJs mask parser (PtrajMask::parseMaskString)
  * which returns a char array of size P->natom, where selected atoms are 
  * denoted by T. Based on this create an array of selected atom numbers. 
  * maskChar is used to determine whether atoms denoted by 'T' or 'F' will
  * be selected (the latter is the case e.g. with stripped atoms). 
  */
int AtomMask::SetupMask(AmberParm *Pin, double *Xin, int debug) {
  char *mask;

  if (Pin==NULL) {
    mprinterr("    Error: AtomMask::SetupMask: (%s) Topology is NULL.\n", maskString.c_str());
    return 1;
  }
  if (Postfix.empty()) {
    mprinterr("    Error: AtomMask::SetupMask: Postfix is NULL.\n");
    return 1;
  }

  if (debug>1) mprintf("\tAtomMask: Parsing postfix [%s]\n",Postfix.c_str());

  // Allocate character atom mask
  mask = parseMaskString((char*)Postfix.c_str(), Pin->natom, Pin->nres, Pin->names, Pin->resnames,
                         Pin->resnums, Xin, Pin->types, debug);
  if (mask==NULL) {
    mprinterr("    Error: Could not set up mask %s for topology %s\n",
              maskString.c_str(), Pin->parmName);
    return 1;
  }

  // Wipe out previous char mask if allocated
  CharMask.clear();

  // Set up an integer list of the selected atoms. 
  // NOTE: For large selections this will use 4x the memory of the char atom
  //       mask. Could check to see which will be bigger.
  Nselected=0;
  Selected.clear();
  for (int atom=0; atom<Pin->natom; atom++) {
    if (mask[atom]==maskChar) {
      Selected.push_back( atom );
      Nselected++;
    }
  }

  if (debug>0) {
    if (maskChar=='F')
      mprintf("          Inverse of Mask %s corresponds to %i atoms.\n",
              maskString.c_str(), Pin->natom - Nselected);
    else
      mprintf("          Mask %s corresponds to %i atoms.\n",maskString.c_str(),
              Nselected);
  }

  // Free the character mask, no longer needed.
  delete[] mask;
  
  return 0;
}

// AtomMask::SetupCharMask()
/** For cases where we need to know both atoms in and out of mask use
  * the output of parseMaskString (1 char for each atom, T for selected,
  * F for not.). 
  */
int AtomMask::SetupCharMask(AmberParm *Pin, double *Xin, int debug) {
  char *mask;
  if (Pin==NULL) {
    mprinterr("    Error: AtomMask::SetupCharMask: (%s) Topology is NULL.\n", maskString.c_str());
    return 1;
  }
  if (Postfix.empty()) {
    mprinterr("    Error: AtomMask::SetupCharMask: Postfix is NULL.\n");
    return 1;
  }

  // Wipe out previous Selected mask if allocated
  Selected.clear();

  // Allocate atom mask - free mask if already allocated
  CharMask.clear();
  mask = parseMaskString((char*)Postfix.c_str(), Pin->natom, Pin->nres, Pin->names, Pin->resnames,
                         Pin->resnums, Xin, Pin->types, debug);
  if (mask==NULL) {
    mprinterr("    Error: Could not set up mask %s for topology %s\n",
              maskString.c_str(), Pin->parmName);
    return 1;
  }
  //CharMask.assign( mask );
  CharMask.reserve( Pin->natom );
  for (char *ptr = mask; *ptr != '\0'; ptr++)
    CharMask.push_back( *ptr );
  delete[] mask;
  // Determine number of selected atoms
  Nselected=0;
  for (int atom=0; atom<Pin->natom; atom++)
    if (CharMask[atom]==maskChar) Nselected++;
  return 0;
}

// AtomMask::AtomInCharMask()
bool AtomMask::AtomInCharMask(int atom) {
  if (CharMask.empty()) return false;
  if (atom < 0) return false;
  if (atom >= (int)CharMask.size()) return false;
  if (CharMask[atom]==maskChar) return true;
  return false;
}

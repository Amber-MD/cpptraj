#include <algorithm> // sort, unique
#include "AtomMask.h"
#include "CpptrajStdio.h"
#include "PtrajMask.h"

// CONSTRUCTOR
AtomMask::AtomMask() {
  Nselected=0;
  Natom = 0;
  maskChar = 'T';
}

// DESTRUCTOR
AtomMask::~AtomMask() {
}

// COPY CONSTRUCTOR
AtomMask::AtomMask(const AtomMask &rhs) {
  Nselected=rhs.Nselected;
  Natom=rhs.Natom;
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
  Natom=rhs.Natom;
  maskChar=rhs.maskChar;
  maskString=rhs.maskString;
  Selected=rhs.Selected;
  CharMask=rhs.CharMask;
  Postfix=rhs.Postfix;
  // Return *this
  return *this;
}

// AtomMask::PostfixExpression()
char *AtomMask::PostfixExpression() { 
  if (Postfix.empty())
    return NULL;
  return (char*)Postfix.c_str(); 
}

// AtomMask::ResetMask()
void AtomMask::ResetMask() {
  Nselected = 0;
  Natom=0;
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

// AtomMask::NumAtomsInCommon()
/** Given an atom mask, determine how many selected atoms this mask
  * has in common with it.
  */
int AtomMask::NumAtomsInCommon(AtomMask &maskIn) {
  std::vector<int> intersect;
  std::vector<int>::iterator intersect_end;

  if (Selected.empty() || maskIn.Selected.empty()) return 0;
  // Max size of the intersection is the min size of either array
  intersect.resize( Selected.size() );
  // Create copies of arrays so they can be sorted
  std::vector<int> selected_1 = Selected;
  std::vector<int> selected_2 = maskIn.Selected;
  // Sort the arrays
  sort(selected_1.begin(), selected_1.end());
  sort(selected_2.begin(), selected_2.end());
  // Set intersect to the intersection of selected_1 and selectd_2
  intersect_end = set_intersection(selected_1.begin(), selected_1.end(),
                                   selected_2.begin(), selected_2.end(),
                                   intersect.begin());
  // DEBUG:
  //mprintf("DBG:\tIntersection of [%s] and [%s] is:",maskString.c_str(),maskIn.maskString.c_str());
  //for ( std::vector<int>::iterator atom = intersect.begin();
  //                                 atom != intersect_end;
  //                                 atom++)
  //  mprintf(" %i",*atom);
  //mprintf("\n");
  return int(intersect_end - intersect.begin());
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
      Nselected = (int) Selected.size();
      return;
    }
  }

  // Add atom to mask
  Selected.push_back(atomIn);
  Nselected = (int) Selected.size();
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
  Nselected = (int) Selected.size();
}

// AtomMask::AddAtomRange()
/** Add atoms in range from minAtom up to but not including maxAtom to 
  * Selected array. The resulting array is sorted and duplicates are removed.
  */
void AtomMask::AddAtomRange(int minAtom, int maxAtom) {
  //mprintf("DEBUG:\t\tAdding atoms %i to %i\n",minAtom,maxAtom);
  if (minAtom >= maxAtom) return;
  for (int atom = minAtom; atom < maxAtom; atom++)
    Selected.push_back( atom );
  // Sort Selected
  sort( Selected.begin(), Selected.end() );
  // Remove duplicates
  std::vector<int>::iterator atomit = unique( Selected.begin(), Selected.end() );
  Selected.resize( atomit - Selected.begin() );
  //mprintf("\t\t[");
  //for (std::vector<int>::iterator da = Selected.begin(); da != Selected.end(); da++)
  //  mprintf(" %i",*da);
  //mprintf("]\n");
  Nselected = (int) Selected.size();
}

// AtomMask::PrintMaskAtoms()
void AtomMask::PrintMaskAtoms(const char *header) {
  mprintf("%s=",header);
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
  mprintf("\n");
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
/** Set up an atom mask containing selected atom numbers given a char
  * array of size natom with T for selected atoms and F for unselected
  * atoms. The actual mask parser is called from AmberParm. 
  * maskChar is used to determine whether atoms denoted by 'T' or 'F' will
  * be selected (the latter is the case e.g. with stripped atoms). 
  */
void AtomMask::SetupMask(char *charmask,int natom,int debug) {
  // Wipe out previous char mask if allocated
  CharMask.clear();

  // Set up an integer list of the selected atoms. 
  // NOTE: For large selections this will use 4x the memory of the char atom
  //       mask. Could check to see which will be bigger.
  Nselected=0;
  Natom = natom;
  Selected.clear();
  for (int atom=0; atom < natom; atom++) {
    if (charmask[atom]==maskChar) {
      Selected.push_back( atom );
      ++Nselected;
    }
  }

  if (debug>0) {
    if (maskChar=='F')
      mprintf("          Inverse of Mask %s corresponds to %i atoms.\n",
              maskString.c_str(), natom - Nselected);
    else
      mprintf("          Mask %s corresponds to %i atoms.\n",maskString.c_str(),
              Nselected);
  }
}

// AtomMask::SetupCharMask()
/** Given an input char mask of size natom, set up a corresponding char mask.
  * Useful for cases where we need to know both atoms in and out of mask.
  */
void AtomMask::SetupCharMask(char *charmask, int natom, int debug) {
  // Wipe out previous Selected mask if allocated
  Selected.clear();

  // Allocate atom mask - free mask if already allocated
  CharMask.clear();

  Nselected=0;
  Natom = natom;
  CharMask.reserve( natom );
  for (int i = 0; i < natom; i++) {
    CharMask.push_back( charmask[i] );
    // Determine number of selected atoms
    if (charmask[i] == maskChar) ++Nselected;
  }
}

// AtomMask::AtomInCharMask()
bool AtomMask::AtomInCharMask(int atom) {
  if (CharMask.empty()) return false;
  if (atom < 0) return false;
  if (atom >= (int)CharMask.size()) return false;
  if (CharMask[atom]==maskChar) return true;
  return false;
}

// AtomMask::ConvertMaskType()
/** If the mask is an integer mask, convert it to a char mask and 
  * vice versa.
  */
int AtomMask::ConvertMaskType() {
  // Nselected remains the same.
  // Integer to Char
  if (!Selected.empty()) {
    // If Natom is 0 (which can happen if mask was set up with AddAtomX
    // functions) this cant work.
    if (Natom==0) {
      mprinterr("Error: AtomMask::ConvertMaskType(): Natom for integer mask is 0.\n");
      return 1;
    }
    CharMask.assign( Natom, 'F' );
    for (std::vector<int>::iterator maskatom = Selected.begin();
                                    maskatom != Selected.end();
                                    maskatom++)
    {
      CharMask[*maskatom]='T';
    }
    Selected.clear();

  // Char to Integer
  } else if (!CharMask.empty()) {
    Selected.reserve( Nselected );
    for (int atom = 0; atom < Natom; atom++) {
      if (CharMask[atom]=='T')
        Selected.push_back( atom );
    }
    CharMask.clear();

  } else {
    mprinterr("Error: AtomMask::ConvertMaskType(): Mask has not been set up.\n");
    return 1;
  }
  return 0;
}


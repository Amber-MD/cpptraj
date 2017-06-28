#include <algorithm> // sort, unique
#include "AtomMask.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
AtomMask::AtomMask(int beginAtom, int endAtom) : Natom_(0), maskChar_(SelectedChar_)
{
  AddAtomRange(beginAtom, endAtom);
}

// CONSTRUCTOR
AtomMask::AtomMask(int atomNum) : Selected_(1, atomNum), Natom_(1), maskChar_(SelectedChar_) {}

// AtomMask::ResetMask()
void AtomMask::ResetMask() {
  Natom_ = 0;
  Selected_.clear();
  ClearTokens();
}

/** Flip the current character used to select atoms. Useful when you want 
  * the mask to select the inverse of the given expression, like in 'strip'.
  */
void AtomMask::InvertMaskExpression() {
  if (maskChar_ == SelectedChar_)
    maskChar_ = UnselectedChar_;
  else
    maskChar_ = SelectedChar_;
}

// AtomMask::InvertMask()
void AtomMask::InvertMask() {
  // Invert the integer mask.
  std::vector<int> invert;
  invert.reserve( Natom_ - (int)Selected_.size() );
  const_iterator selected_atom = Selected_.begin();
  for (int idx = 0; idx < Natom_; idx++) {
    if (selected_atom == Selected_.end() || idx != *selected_atom) {
      // Atom was not selected or no more selected atoms; add.
      invert.push_back( idx );
    } else {
      // Atom was selected; ignore and advance to next selected atom.
      ++selected_atom;
    }
  }
  Selected_ = invert;
}

// AtomMask::NumAtomsInCommon()
/** Given an atom mask, determine how many selected atoms this mask
  * has in common with it.
  */
int AtomMask::NumAtomsInCommon(AtomMask const& maskIn) {
  std::vector<int> intersect;
  std::vector<int>::iterator intersect_end;

  if (Selected_.empty() || maskIn.Selected_.empty()) return 0;
  // Max size of the intersection is the min size of either array
  intersect.resize( Selected_.size() );
  // Create copies of arrays so they can be sorted
  std::vector<int> selected_1 = Selected_;
  std::vector<int> selected_2 = maskIn.Selected_;
  // Sort the arrays
  std::sort(selected_1.begin(), selected_1.end());
  std::sort(selected_2.begin(), selected_2.end());
  // Set intersect to the intersection of selected_1 and selectd_2
  intersect_end = std::set_intersection(selected_1.begin(), selected_1.end(),
                                   selected_2.begin(), selected_2.end(),
                                   intersect.begin());
  // DEBUG:
  //mprintf("DBG:\tIntersection of [%s] and [%s] is:",maskString_.c_str(),maskIn.maskString_.c_str());
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
  for (std::vector<int>::iterator atom = Selected_.begin(); atom != Selected_.end(); atom++) {
    if ( *atom == atomIn) return;
    if ( *atom > atomIn) {
      // Insert at the current position, which is the first atom # > atomIn
      Selected_.insert(atom, atomIn);
      return;
    }
  }
  // Add atom to mask
  Selected_.push_back(atomIn);
}

// AtomMask::AddAtoms()
/** Given an array, add the atom numbers in array to the Selected_ array.
  * The resulting array is sorted and any duplicates are removed.
  */
void AtomMask::AddAtoms(std::vector<int> const& atomsIn) {
  std::vector<int>::const_iterator atom;
  // Make room for atomsIn in Selected_
  //Selected_.reserve( Selected_.size() + atomsIn.size() );
  // Put every atom in atomsIn in Selected_ array
  for (atom = atomsIn.begin(); atom != atomsIn.end(); atom++) 
    Selected_.push_back( *atom ); 
  // Sort Selected_
  std::sort( Selected_.begin(), Selected_.end() );
  // Remove duplicates
  atom = unique( Selected_.begin(), Selected_.end() );
  Selected_.resize( atom - Selected_.begin() );
}

// AtomMask::AddAtomRange()
/** Add atoms in range from minAtom up to but not including maxAtom to 
  * Selected_ array. The resulting array is sorted and duplicates are removed.
  */
void AtomMask::AddAtomRange(int minAtom, int maxAtom) {
  //mprintf("DEBUG:\t\tAdding atoms %i to %i\n",minAtom,maxAtom);
  if (minAtom >= maxAtom) return;
  for (int atom = minAtom; atom < maxAtom; atom++)
    Selected_.push_back( atom );
  // Sort Selected_
  std::sort( Selected_.begin(), Selected_.end() );
  // Remove duplicates
  std::vector<int>::iterator atomit = unique( Selected_.begin(), Selected_.end() );
  Selected_.resize( atomit - Selected_.begin() );
  //mprintf("\t\t[");
  //for (std::vector<int>::iterator da = Selected_.begin(); da != Selected_.end(); da++)
  //  mprintf(" %i",*da);
  //mprintf("]\n");
}

// AtomMask::AddMaskAtPosition()
/** Add the atoms in Mask to this mask starting at the specified positon.
  * Currently only used by Closest for modifying the original stripped
  * mask with the calculated closest waters.
  */
void AtomMask::AddMaskAtPosition(AtomMask const& maskIn, int idx) {
  // NOTE: NO BOUNDS CHECK!
  for (const_iterator atom = maskIn.begin(); atom != maskIn.end(); atom++) 
    Selected_[idx++] = *atom;
}

// AtomMask::PrintMaskAtoms()
void AtomMask::PrintMaskAtoms(const char *header) const {
  mprintf("%s=",header);
  if (!Selected_.empty()) {
    for (std::vector<int>::const_iterator atom = Selected_.begin(); 
                                          atom != Selected_.end(); ++atom)
      mprintf(" %i",*atom + 1);
  } else 
    mprintf("No atoms selected.");
  mprintf("\n");
}

/** Set up an atom mask containing selected atom numbers given a char
  * array of size natom with T for selected atoms and F for unselected
  * atoms. The actual mask parser is called from Topology. 
  * maskChar_ is used to determine whether atoms denoted by 'T' or 'F' will
  * be selected (the latter is the case e.g. with stripped atoms). 
  */
int AtomMask::SetupMask(AtomArrayT const& atoms, ResArrayT const& residues, const double* XYZ)
{
  // Set up an integer list of the selected atoms. 
  // NOTE: For large selections this will use 4x the memory of the char atom
  //       mask. Could check to see which will be bigger.
  Natom_ = (int)atoms.size();
  Selected_.clear();
  char* charmask = ParseMask(atoms, residues, XYZ);
  if (charmask == 0) return 1;
  for (int atom = 0; atom != Natom_; atom++) {
    if (charmask[atom] == maskChar_)
      Selected_.push_back( atom );
  }
  delete[] charmask;
  return 0;
}

// AtomMask::ConvertToCharMask()
/** Can be used to set up a CharMask, e.g.
  * CharMask mask( AtomMask.ConvertToCharMask(), AtomMask.Nselected() )
  */
std::vector<char> AtomMask::ConvertToCharMask() const {
  // If Natom is empty this will not work.
  if (Natom_ < 1) {
    mprinterr("Internal Error: Cannot convert AtomMask to CharMask, Natom < 1.\n");
    return std::vector<char>();
  }
  std::vector<char> CharMask(Natom_, UnselectedChar_);
  if (!Selected_.empty()) {
    for (std::vector<int>::const_iterator maskatom = Selected_.begin();
                                          maskatom != Selected_.end(); ++maskatom)
      CharMask[*maskatom] = SelectedChar_;
  }
  return CharMask;
}

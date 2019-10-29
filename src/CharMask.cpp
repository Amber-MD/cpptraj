#include "CharMask.h"
#include "CpptrajStdio.h" // PrintMaskAtoms

/** Use to initialize the mask without a mask expression. */
void CharMask::InitCharMask(int natoms, bool initSelected) {
  if (initSelected)
    CharMask_.assign(natoms, SelectedChar_);
  else
    CharMask_.assign(natoms, UnselectedChar_);
}

/** Given atom and residue info and coordinates, setup character mask
  * based on current mask tokens.
  */
int CharMask::SetupMask(AtomArrayT const& atoms, ResArrayT const& residues,
                        MolArrayT const& molecules, const double* XYZ)
{
  CharMask_.clear();
  nselected_ = 0;
  CharMask_.reserve( atoms.size() );
  char* charmask = ParseMask(atoms, residues, molecules, XYZ);
  if (charmask == 0) return 1;
  for (unsigned int i = 0; i != atoms.size(); i++) {
    CharMask_.push_back( charmask[i] );
    if (charmask[i] == SelectedChar_) ++nselected_;
  }
  delete[] charmask;
  return 0;
}

// CharMask::PrintMaskAtoms()
void CharMask::PrintMaskAtoms(const char *header) const {
  mprintf("%s=",header);
  if (!CharMask_.empty()) {
    for (unsigned int atomnum = 0; atomnum != CharMask_.size(); ++atomnum)
      if (CharMask_[atomnum] == SelectedChar_) mprintf(" %i", atomnum+1);
  } else
    mprintf(" No atoms selected.");
  mprintf("\n");
}

// CharMask::ResetMask()
void CharMask::ResetMask() {
  CharMask_.clear();
  nselected_ = 0;
  ClearTokens();
}

// CharMask::ClearSelected()
void CharMask::ClearSelected() {
  CharMask_.assign(CharMask_.size(), UnselectedChar_);
  nselected_ = 0;
}

// CharMask::InvertMask()
void CharMask::InvertMask() {
  for (std::vector<char>::iterator atchar = CharMask_.begin();
                                   atchar != CharMask_.end(); ++atchar)
    if ( *atchar == SelectedChar_ )
      *atchar = UnselectedChar_;
     else
      *atchar = SelectedChar_;
  nselected_ = (int)CharMask_.size() - nselected_;
}

// CharMask::AtomInCharMask()
bool CharMask::AtomInCharMask(int atom) const {
  if (CharMask_.empty()) return false;
  if (atom < 0) return false;
  if (atom >= (int)CharMask_.size()) return false;
  if (CharMask_[atom] == SelectedChar_) return true;
  return false;
}

// CharMask::AtomsInCharMask()
bool CharMask::AtomsInCharMask(int startatom, int endatom) const {
  if (CharMask_.empty()) return false;
  if (startatom > endatom) return false;
  if (startatom < 0) return false;
  if (endatom > (int)CharMask_.size()) return false;
  for (int idx = startatom; idx < endatom; ++idx)
    if (CharMask_[idx] == SelectedChar_) return true;
  return false;
}

/** This routine can be used to convert a CharMask to an AtomMask, e.g.
  * AtomMask mask( CharMask.ConvertToIntMask(), CharMask.Natom() )
  */
std::vector<int> CharMask::ConvertToIntMask() const {
  std::vector<int> Selected;
  if (CharMask_.empty()) return Selected;
  Selected.reserve( nselected_ );
  for (int atom = 0; atom != (int)CharMask_.size(); atom++) {
    if (CharMask_[atom] == SelectedChar_)
      Selected.push_back( atom );
  }
  return Selected;
}

#include "CharMask.h"
#include "CpptrajStdio.h" // PrintMaskAtoms

int CharMask::SetupMask(AtomArrayT const& atoms, ResArrayT const& residues, const double* XYZ)
{
  CharMask_.clear();
  nselected_ = 0;
  CharMask_.reserve( atoms.size() );
  char* charmask = ParseMask(atoms, residues, XYZ);
  if (charmask == 0) return 1;
  for (unsigned int i = 0; i != atoms.size(); i++) {
    CharMask_.push_back( charmask[i] );
    if (charmask[i] == maskChar_) ++nselected_;
  }
  delete[] charmask;
  return 0;
}

void CharMask::PrintMaskAtoms(const char *header) const {
  mprintf("%s=",header);
  if (!CharMask_.empty()) {
    for (unsigned int atomnum = 0; atomnum != CharMask_.size(); ++atomnum)
      if (CharMask_[atomnum] == maskChar_) mprintf(" %i", atomnum+1);
  } else
    mprintf(" No atoms selected.");
  mprintf("\n");
}

void CharMask::ResetMask() {
  CharMask_.clear();
  nselected_ = 0;
  ClearTokens();
}

void CharMask::ClearSelected() {
  for (std::vector<char>::iterator m = CharMask_.begin(); m != CharMask_.end(); ++m)
    *m = unselectedChar_;
  nselected_ = 0;
}

void CharMask::InvertMask() {
  for (std::vector<char>::iterator atchar = CharMask_.begin();
                                   atchar != CharMask_.end(); ++atchar)
    if ( *atchar == maskChar_ )
      *atchar = unselectedChar_;
     else
      *atchar = maskChar_;
  nselected_ = (int)CharMask_.size() - nselected_;
}

bool CharMask::AtomInCharMask(int atom) const {
  if (CharMask_.empty()) return false;
  if (atom < 0) return false;
  if (atom >= (int)CharMask_.size()) return false;
  if (CharMask_[atom] == maskChar_) return true;
  return false;
}

bool CharMask::AtomsInCharMask(int startatom, int endatom) const {
  if (CharMask_.empty()) return false;
  if (startatom > endatom) return false;
  if (startatom < 0) return false;
  if (endatom > (int)CharMask_.size()) return false;
  for (int idx = startatom; idx < endatom; ++idx)
    if (CharMask_[idx] == maskChar_) return true;
  return false;
}

std::vector<int> CharMask::ConvertToIntMask() const {
  std::vector<int> Selected;
  if (CharMask_.empty()) return Selected;
  Selected.reserve( nselected_ );
  for (int atom = 0; atom != (int)CharMask_.size(); atom++) {
    if (CharMask_[atom] == maskChar_)
      Selected.push_back( atom );
  }
  return Selected;
}

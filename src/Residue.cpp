#include "Residue.h"

Residue::Residue() :
  original_resnum_(0),
  firstAtom_(0),
  resname_("")
{}

Residue::Residue(int resnum, NameType resname) :
  original_resnum_(resnum), firstAtom_(0), resname_(resname)
{}

Residue::Residue(NameType resname, int firstAtomIn) :
  original_resnum_(0), firstAtom_(firstAtomIn), resname_(resname) 
{}

Residue::Residue(int firstAtomIn) :
  original_resnum_(-1),
  firstAtom_(firstAtomIn),
  resname_("")
{} 

bool Residue::operator==(const Residue &rhs) {
  return (original_resnum_ == rhs.original_resnum_);
}

bool Residue::operator!=(const Residue &rhs) {
  return (original_resnum_ != rhs.original_resnum_);
}

void Residue::SetFirstAtom(int firstAtomIn) {
  firstAtom_ = firstAtomIn;
}

bool Residue::NameIsSolvent() {
  if (resname_ == "WAT " ||
      resname_ == "HOH " ||
      resname_ == "TIP3"
     )
    return true;
  return false;
}


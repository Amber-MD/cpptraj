#include "Molecule.h"

Molecule::Molecule(): 
  firstResNum_(0),
  beginAtom_(0),
  endAtom_(0),
  isSolvent_(false)
{ }

Molecule::Molecule(int firstRes, int begin, int end) :
  firstResNum_(firstRes),
  beginAtom_(begin),
  endAtom_(end),
  isSolvent_(false)
{ }

void Molecule::SetFirst(int begin, int firstRes) {
  firstResNum_ = firstRes;
  beginAtom_ = begin;
}

void Molecule::SetLast(int end) {
  endAtom_ = end;
}

void Molecule::SetSolvent() {
  isSolvent_ = true;
}


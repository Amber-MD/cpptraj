#include "Molecule.h"

Molecule::Molecule(): 
  firstResNum_(0),
  beginAtom_(0),
  endAtom_(0)
{ }

Molecule::Molecule(int firstRes, int begin, int end) :
  firstResNum_(firstRes),
  beginAtom_(begin),
  endAtom_(end)
{ }

void Molecule::SetFirst(int begin, int firstRes) {
  firstResNum_ = firstRes;
  beginAtom_ = begin;
}

void Molecule::SetLast(int end) {
  endAtom_ = end;
}


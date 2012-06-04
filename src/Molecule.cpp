#include "Molecule.h"

Molecule::Molecule(): 
  beginAtom_(0),
  endAtom_(0),
  isSolvent_(false)
{ }

Molecule::Molecule(int begin, int end) :
  beginAtom_(begin),
  endAtom_(end),
  isSolvent_(false)
{ }


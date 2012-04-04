#include "MapAtom.h"

const char MapAtom::AtomicElementChar[Atom::NUMELEMENTS] = { 0,
  'H',  'B',  'C',  'N',  'O', 'F',  'P',  'S',  'X',  'Y',  'f',  'c',
  'I',  'M',  'U',  'L',  'K', 'R',  'E',  'Z',  'n'
};

MapAtom::MapAtom() :
  isChiral_(false), 
  isMapped_(false), 
  complete_(false), 
  Nduplicated_(0),
  name_(0) 
{}

MapAtom::MapAtom(const MapAtom& rhs) : Atom(rhs),
   isChiral_(rhs.isChiral_), 
   isMapped_(rhs.isMapped_), 
   complete_(rhs.complete_),
   atomID_(rhs.atomID_), 
   unique_(rhs.unique_), 
   Nduplicated_(rhs.Nduplicated_),
   name_(rhs.name_) 
{}

MapAtom::MapAtom(const Atom& rhs) : Atom(rhs),
  isChiral_(false),
  isMapped_(false),
  complete_(false),
  Nduplicated_(0),
  name_(AtomicElementChar[Element()])
{}

MapAtom& MapAtom::operator=(const MapAtom& rhs) {
  if (&rhs == this) return *this;
  Atom::operator=(rhs);
  isChiral_ = rhs.isChiral_;
  isMapped_ = rhs.isMapped_;
  complete_ = rhs.complete_;
  atomID_   = rhs.atomID_;
  unique_   = rhs.unique_;
  Nduplicated_ = rhs.Nduplicated_;
  name_ = rhs.name_;
  return *this;
}


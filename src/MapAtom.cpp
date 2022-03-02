#include "MapAtom.h"

/** Try to assign more common elements letters that are like
  * their names; for less common elements just use whatever
  * characters are handy. 
  */
const char MapAtom::AtomicElementChar_[Atom::NUMELEMENTS_] = { 0,
  'H',  'B',  'C',  'N',  'O', 'F',  
  'P',  'S',  'X',  'Y',  'f', 'c',
  'I',  'M',  'U',  'L',  'K', 'R',  
  'E',  'Z',  'n',  'A',  'r', 'a',
  's',  'G',  't',  'b',  '{', 'i',
  'h',  'o',  'D',  '}',  '|', '~',
  '!',  '"',  '#',  '$',  '%', '&',
  '(',  ')',  '*',  '+',  ',', '-',
  '.',  '/',  '0',  '1',  '2', '3',
  '4',  '5',  '6',  '7',  '8', '9',
  ':',  ';',  '<',  '=',  '>', '?',
  '@',  '[', '\\',  ']',  '^', '_',
  '\'',  '`', 'm',
  'x'
};

/// CONSTRUCTOR
MapAtom::MapAtom() :
  isChiral_(false),
  boundToChiral_(false), 
  isMapped_(false), 
  complete_(false), 
  Nduplicated_(0),
  name_(0) 
{
  std::fill( xyz_, xyz_ + 3, 0.0 );
}

// COPY CONSTRUCTOR
MapAtom::MapAtom(const MapAtom& rhs) : Atom(rhs),
   isChiral_(rhs.isChiral_),
   boundToChiral_(rhs.boundToChiral_), 
   isMapped_(rhs.isMapped_), 
   complete_(rhs.complete_),
   atomID_(rhs.atomID_), 
   unique_(rhs.unique_), 
   Nduplicated_(rhs.Nduplicated_),
   name_(rhs.name_) 
{
  std::copy( rhs.xyz_, rhs.xyz_ + 3, xyz_ );
}

// COPY CONSTRUCTOR
/// Copy base atom to this MapAtom
MapAtom::MapAtom(const Atom& rhs, const double* xyzIn) : Atom(rhs),
  isChiral_(false),
  boundToChiral_(false),
  isMapped_(false),
  complete_(false),
  Nduplicated_(0),
  name_(AtomicElementChar_[Element()])
{
  std::copy( xyzIn, xyzIn + 3, xyz_ );
}

// Assignment
MapAtom& MapAtom::operator=(const MapAtom& rhs) {
  if (&rhs == this) return *this;
  Atom::operator=(rhs);
  isChiral_ = rhs.isChiral_;
  boundToChiral_ = rhs.boundToChiral_;
  isMapped_ = rhs.isMapped_;
  complete_ = rhs.complete_;
  atomID_   = rhs.atomID_;
  unique_   = rhs.unique_;
  Nduplicated_ = rhs.Nduplicated_;
  name_ = rhs.name_;
  std::copy( rhs.xyz_, rhs.xyz_ + 3, xyz_ );
  return *this;
}

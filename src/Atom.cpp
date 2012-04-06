#include <algorithm>
#include "Atom.h"
#include "CpptrajStdio.h"

const int Atom::AtomicElementNum[NUMELEMENTS] = { 0,
 1,  3,  6,  7, 8,  9,  15, 16, 17, 35, 26, 20, 
 53, 12, 29, 3, 19, 37, 55, 30, 11
};

/// Atom names corresponding to AtomicElementType.
// 2 chars + NULL.
const char Atom::AtomicElementName[NUMELEMENTS][3] = { "??",
  "H",  "B",  "C",  "N", "O",  "F",  "P",  "S", "Cl", "Br", "Fe", "Ca",
  "I", "Mg", "Cu", "Li", "K", "Rb", "Cs", "Zn", "Na"
};

// CONSTRUCTOR
Atom::Atom() : 
  charge_(0),
  mass_(1),
  gb_radius_(0),
  gb_screen_(0),
  aname_(""),
  atype_(""),
  atype_index_(0),
  element_(UNKNOWN_ELEMENT),
  resnum_(0),
  mol_(0)
{ }

// CONSTRUCTOR
/// Take atom name and coordinates. Attempt to determine element from name.
//Atom::Atom(NameType aname, double (&XYZ)[3]) :
Atom::Atom(NameType aname) :
  charge_(0),
  mass_(1),
  gb_radius_(0),
  gb_screen_(0),
  aname_(aname),
  atype_(""),
  atype_index_(0),
  element_(UNKNOWN_ELEMENT),
  resnum_(0),
  mol_(0)
{
  SetElementFromName();
}

// CONSTRUCTOR
//Atom::Atom( NameType aname, double (&XYZ)[3], NameType atype, double q ) :
Atom::Atom( NameType aname, NameType atype, double q ) :
  charge_(q),
  mass_(1),
  gb_radius_(0),
  gb_screen_(0),
  aname_(aname),
  atype_(atype),
  atype_index_(0),
  element_(UNKNOWN_ELEMENT),
  resnum_(0),
  mol_(0)
{
  SetElementFromName();
}

// CONSTRUCTOR
Atom::Atom( NameType name, double charge, int atomicnum, double mass, int atidx,
            NameType type, double rad, double screen, int resnum ) :
  charge_(charge),
  mass_(mass),
  gb_radius_(rad),
  gb_screen_(screen),
  aname_(name),
  atype_(type),
  atype_index_(atidx),
  element_(UNKNOWN_ELEMENT),
  resnum_(resnum),
  mol_(0)
{
  // Determine element from atomic number
  if (atomicnum>0) {
    for (int i = 1; i < 22; i++)
      if (AtomicElementNum[i] == atomicnum)
        element_ = (AtomicElementType) i;
  }
  // If element still unknown attempt to determine from name
  if (element_ == UNKNOWN_ELEMENT)
    SetElementFromName();
}

// COPY CONSTRUCTOR
Atom::Atom(const Atom &rhs) :
  charge_(rhs.charge_),
  mass_(rhs.mass_),
  gb_radius_(rhs.gb_radius_),
  gb_screen_(rhs.gb_screen_),
  aname_(rhs.aname_),
  atype_(rhs.atype_),
  atype_index_(rhs.atype_index_),
  element_(rhs.element_),
  resnum_(rhs.resnum_),
  mol_(rhs.mol_),
  bonds_(rhs.bonds_)
{ }

// SWAP
void Atom::swap(Atom &first, Atom &second) {
  using std::swap;
  swap(first.charge_, second.charge_);
  swap(first.mass_, second.mass_);
  swap(first.gb_radius_, second.gb_radius_);
  swap(first.gb_screen_, second.gb_screen_);
  swap(first.aname_, second.aname_);
  swap(first.atype_, second.atype_);
  swap(first.atype_index_, second.atype_index_);
  swap(first.element_, second.element_);
  swap(first.resnum_, second.resnum_);
  swap(first.mol_, second.mol_);
  swap(first.bonds_, second.bonds_);
}

// ASSIGNMENT via copy/swap idiom
Atom &Atom::operator=(Atom other) {
  swap(*this, other);
  return *this;
}

// Atom::Info()
void Atom::Info() {
  mprintf(" [%s]",*aname_);
  mprintf(" Res %i:",resnum_+1);
  mprintf(" Mol %i", mol_+1);
  mprintf(" Type=[%s]",*atype_);
  mprintf(" Charge=%lf",charge_);
  mprintf(" Mass=%lf",mass_);
  mprintf("\n");
}

// Atom::SetName()
void Atom::SetName(NameType nameIn) {
  aname_ = nameIn;
}

// Atom::SetResNum()
void Atom::SetResNum(int resnumIn) {
  resnum_ = resnumIn;
}

// Atom::SetMol()
void Atom::SetMol(int molIn) {
  mol_ = molIn;
}

// Atom::NoMol()
bool Atom::NoMol() {
  return ( mol_ < 0 );
}

// Atom::AddBond()
void Atom::AddBond(int idxIn) {
  bonds_.push_back( idxIn );
}

// Atom::ClearBonds()
void Atom::ClearBonds() {
  bonds_.clear();
}

void Atom::SortBonds() {
  sort( bonds_.begin(), bonds_.end() );
}

// Atom::AddExcluded()
/*void Atom::AddExcluded(int idxIn) {
  excluded_.insert( idxIn );
}*/
void Atom::AddExclusionList(std::set<int>& elist) {
  excluded_.clear();
  for (std::set<int>::iterator ei = elist.begin(); ei != elist.end(); ei++)
    excluded_.push_back( *ei );
}

// Atom::ClearExcluded()
/*void Atom::ClearExcluded() {
  excluded_.clear();
}*/

// Atom::SetElementFromName()
/** If not already known, try to determine atomic element from atom name. 
  * Based on Amber standard atom names.
  */
void Atom::SetElementFromName() {
  if (element_!=UNKNOWN_ELEMENT) return;
  char c1 = aname_[0];
  if (c1=='\0') return;
  char c2 = aname_[1];

  switch (c1) {
    case 'H' : element_ = HYDROGEN; break;
    case 'C' :
      if (c2=='l' || c2=='L') element_ = CHLORINE;
      else if (c2=='0') element_ = CALCIUM;
      else if (c2=='s') element_ = CESIUM;
      else if (c2=='U') element_ = COPPER;
      else element_ = CARBON;
      break;
    case 'N' :
      if (c2=='a') element_ = SODIUM;
      else element_ = NITROGEN;
      break;
    case 'O' : element_ = OXYGEN; break;
    case 'F' :
      if (c2=='e' || c2=='E') element_ = IRON;
      else element_ = FLUORINE;
      break;
    case 'B' :
      if (c2=='r' || c2=='R') element_ = BROMINE;
      else element_ = BORON;
      break;
    case 'I' :
      if (c2=='M') element_ = CHLORINE;    // IM, Cl- in old FF94
      else if (c2=='P') element_ = SODIUM; // IP, Na+ in old FF94
      else element_ = IODINE;
      break;
    case 'P' : element_ = PHOSPHORUS; break;
    case 'S' : element_ = SULFUR; break;
    case 'M' :
      if (c2=='G') element_ = MAGNESIUM;
      break;
    case 'Z' :
      if (c2=='n') element_ = ZINC;
      break;
    case 'L' :
      if (c2=='i') element_ = LITHIUM;
      break;
    case 'K' : element_ = POTASSIUM; break;
    case 'R' :
      if (c2=='b') element_ = RUBIDIUM;
      break;
    default:
      mprintf("Warning: Could not determine atomic number from name [%s]\n",*aname_);
  }
}


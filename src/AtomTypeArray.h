#ifndef INC_ATOMTYPEARRAY_H
#define INC_ATOMTYPEARRAY_H
#include <vector>
#include <map>
#include "AtomType.h"
/// Map unique atom type names to atom types
class AtomTypeArray {
  public:
    AtomTypeArray() {}

    bool AddAtomType(NameType const&, AtomType const&);
  private:
    typedef std::vector<AtomType> Tarray;
    typedef std::map<NameType, int> Tmap;
    typedef std::pair<NameType, int> Tpair;

    Tarray types_;
    Tmap nameToIdx_;
};
#endif

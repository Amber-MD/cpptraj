#ifndef INC_ATOMTYPEARRAY_H
#define INC_ATOMTYPEARRAY_H
#include <vector>
#include <map>
#include "AtomType.h"
/// Map unique atom type names to atom types
class AtomTypeArray {
    typedef std::vector<AtomType> Tarray;
    typedef std::map<NameType, int> Tmap;
  public:
    AtomTypeArray() {}

    AtomType const& operator[](int idx) const { return types_[idx]; }

    bool AddAtomType(NameType const&, AtomType const&);
    int CheckForAtomType(NameType const& n);

    AtomType& UpdateType(int idx) { return types_[idx]; }

    typedef Tmap::const_iterator const_iterator;
    const_iterator begin() const { return nameToIdx_.begin(); }
    const_iterator end()   const { return nameToIdx_.end();   }
  private:
    typedef std::pair<NameType, int> Tpair;

    Tarray types_;
    Tmap nameToIdx_;
};
#endif

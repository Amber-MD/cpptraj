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
    AtomTypeArray() : debug_(0) {}

    AtomType const& operator[](int idx) const { return types_[idx]; }

    void SetDebug(int d) { debug_ = d; }
    /// \return true if type name already present, false otherwise.
    bool AddAtomType(NameType const&, AtomType const&);
    /// \return Atom type index of new/existing atom type.
    int CheckForAtomType(NameType const&, AtomType const&);
    int CheckForAtomType(NameType const&);
    /// \return Atom type index of given name, -1 if not present.
    int AtomTypeIndex(NameType const&);

    AtomType& UpdateType(int idx) { return types_[idx]; }
    unsigned int Size() const { return types_.size(); }

    typedef Tmap::const_iterator const_iterator;
    const_iterator begin() const { return nameToIdx_.begin(); }
    const_iterator end()   const { return nameToIdx_.end();   }
  private:
    typedef std::pair<NameType, int> Tpair;

    Tarray types_;
    Tmap nameToIdx_;
    int debug_;
};
#endif

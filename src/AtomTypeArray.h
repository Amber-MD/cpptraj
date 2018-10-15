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
    /// \return const reference to type at specified position
    AtomType const& operator[](int idx) const { return types_[idx]; }
    /// Set the debug level
    void SetDebug(int d) { debug_ = d; }
    /// \return true if type name already present, false otherwise.
    bool AddAtomType(NameType const&, AtomType const&);
    /// \return Atom type index of new/existing atom type.
    int CheckForAtomType(NameType const&, AtomType const&);
    /// \return Atom type index of new/existing atom type. New types will be placeholders.
    int CheckForAtomType(NameType const&);
    /// \return Atom type index of given type name, -1 if not present.
    int AtomTypeIndex(NameType const&);
    /// \return Reference to type at specified position
    AtomType& UpdateType(int idx) { return types_[idx]; }
    /// \return Number of atom types in array.
    unsigned int Size() const { return types_.size(); }
    /// Iterator over type name to type index map
    typedef Tmap::const_iterator const_iterator;
    /// \return iterator to beginning of name/index map.
    const_iterator begin() const { return nameToIdx_.begin(); }
    /// \return iterator to end of name/index map.
    const_iterator end()   const { return nameToIdx_.end();   }
  private:
    typedef std::pair<NameType, int> Tpair;

    Tarray types_;   ///< Array of atom types.
    Tmap nameToIdx_; ///< Map atom type name to atom type index.
    int debug_;      ///< Debug level
};
#endif

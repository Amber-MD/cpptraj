#ifndef INC_EXCLUSIONARRAY_H
#define INC_EXCLUSIONARRAY_H
#include <set>
#include <vector>
// Forward declares
class Atom;
class AtomMask;
/// Array holding what atoms should be excluded from interacting with a given atom.
class ExclusionArray {
  public:
    /// CONSTRUCTOR
    ExclusionArray();

    enum SelfOpt {
      EXCLUDE_SELF = 0, ///< Include i in exclusion list for atom i.
      NO_EXCLUDE_SELF   ///< Do not include self in exclusion list.
    };
    enum ListOpt {
      FULL = 0,        ///< Exclusion list contains all excluded atoms.
      ONLY_GREATER_IDX ///< Exclusion list only contains atoms with larger index
    };

    /// Hold list of excluded atoms.
    typedef std::set<int> ExListType;
    /// Hold array of exclusion lists
    typedef std::vector<ExListType> ExArrayType;

    /// Exclusion array iterator
    typedef ExArrayType::const_iterator const_iterator;
    /// \return const iterator to beginning of excluded atom array
    const_iterator begin() const { return Excluded_.begin(); }
    /// \return const iterator to end of excluded atom array
    const_iterator end()   const { return Excluded_.end();   }

    /// Exclusion list iterator
    //typedef ExListType::const_iterator ExListIt;
    /// \return Exclusion list for atom with given index.
    ExListType const& operator[](int idx) const { return Excluded_[idx]; }

    /// Set up exclusion list for specified atoms and distance
    int SetupExcluded(std::vector<Atom> const&, AtomMask const&, int, SelfOpt, ListOpt);
    /// Set up exclusion list for all atoms with specified distance.
    int SetupExcluded(std::vector<Atom> const&, int, SelfOpt, ListOpt);
  private:
    /// Determine through-bond distance between atoms
    static void AtomDistance(std::vector<Atom> const&, int, int, int, ExListType&, int);
    /// Set up exclusion array for specified atom and distance
    //static void DetermineExcludedAtoms(ExListType&, std::vector<Atom> const&, int, int);

    ExArrayType Excluded_; ///< Hold exclusion list for each atom.
};
#endif

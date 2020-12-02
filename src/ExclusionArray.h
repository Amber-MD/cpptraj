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
    typedef ExListType::const_iterator ExListIt;

    /// Set up exclusion list
    int SetupExcluded(std::vector<Atom> const&, AtomMask const&);
  private:
    ExArrayType Excluded_; ///< Hold exclusion list for each atom.
};
#endif

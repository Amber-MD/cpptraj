#ifndef INC_PARM_DIHEDRALPARMHOLDER_H
#define INC_PARM_DIHEDRALPARMHOLDER_H
#include <vector>
#include <utility> // std::pair
#include "ParmEnum.h"
#include "../ParameterTypes.h"
#include "../TypeNameHolder.h"
namespace Cpptraj {
namespace Parm {
// -----------------------------------------------------------------------------
/// Specialized class for associating atom types with dihedral parameters.
/** Dihedrals are specialized because you may need to do wildcard matching.
  */
class DihedralParmHolder {
    typedef std::pair<TypeNameHolder,DihedralParmArray> Bpair;
    typedef std::vector<Bpair> Bmap;
  public:
    /// CONSTRUCTOR
    DihedralParmHolder() {}
    /// DESTRUCTOR - virtual since inherited
    virtual ~DihedralParmHolder() {}
    /// Clear dihedral parameters
    void clear()              { bpmap_.clear();        }
    /// \return Number of dihedral parameters
    size_t size()       const { return bpmap_.size();  }
    /// \return true if no dihedral parameters
    bool empty()        const { return bpmap_.empty(); }
    /// \return Last parameter to be overwritten from AddParm()
    DihedralParmArray const& PreviousParm() const { return previousParm_; }
    /// Set wildcard character
    void SetWildcard(char wc) { wc_ = NameType(std::string(1, wc)); }
    /// Add array of dihedral parameters with unique multiplicities
    RetType AddParm(TypeNameHolder const& types, DihedralParmArray const& dpa, bool allowUpdate) {
      // Check if parm for these types exist
      Bmap::iterator it0 = bpmap_.begin();
      for (; it0 != bpmap_.end(); ++it0)
      {
        if (it0->first.Match_NoWC( types ))
          break;
      }
      if (it0 == bpmap_.end()) {
        // Brand new dihedral for these types.
        bpmap_.push_back( Bpair(types, dpa) );
      } else {
        if (!allowUpdate) return ERR;
        // Check if sizes are the same.
        bool update = false;
        if (it0->second.size() != dpa.size())
          update = true;
        else {
          // Sizes are the same. See if parameters are the same.
          for (unsigned int i = 0; i != it0->second.size(); i++) {
            if (it0->second[i] < dpa[i] || dpa[i] < it0->second[i]) {
              update = true;
              break;
            }
          }
        }
        if (update) {
          previousParm_ = it0->second;
          it0->second = dpa;
          return UPDATED;
        } else
          return SAME;
      }
      return ADDED;
    }
    /// Const iterator
    typedef typename Bmap::const_iterator const_iterator;
    /// \return const iterator to beginning
    const_iterator begin() const { return bpmap_.begin(); }
    /// \return const iterator to end
    const_iterator end()   const { return bpmap_.end();   }
    /// \return Array of dihedral parameters matching given atom types.
    DihedralParmArray FindParam(TypeNameHolder const& types, bool& found) const {
      found = true;
      for (const_iterator it = begin(); it != end(); ++it)
        if (it->first.Match_NoWC( types )) return it->second;
      if (wc_.len() > 0) {
        for (const_iterator it = begin(); it != end(); ++it)
          if (it->first.Match_WC( types, wc_)) return it->second;
      }
      found = false;
      return DihedralParmArray();
    }
    /// \return iterator to parameter matching the given types.
    const_iterator GetParam(TypeNameHolder const& types) {
      for (const_iterator it = begin(); it != end(); ++it)
        if (it->first.Match_NoWC( types )) return it;
      if (wc_.len() > 0) {
        for (const_iterator it = begin(); it != end(); ++it)
          if (it->first.Match_WC( types, wc_)) return it;
      }
      return end();
    }
    /// \return size in memory in bytes
    size_t DataSize() const {
      if (bpmap_.empty()) return 0;
      const_iterator elt0 = begin();
      // Assume all TypeNameHolders are the same size
      return (bpmap_.size() * elt0->first.DataSize()) +
             (bpmap_.size() * sizeof(DihedralParmArray)) +
             sizeof(Bmap);
    }
  protected:
    NameType wc_; ///< Wildcard character
  private:
    Bmap bpmap_;
    DihedralParmArray previousParm_; ///< When parameter is updated, store previous value.
};
} // END namespace Parm
} // END namespace Cpptraj
#endif

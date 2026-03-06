#ifndef INC_PARM_PARMHOLDER_H
#define INC_PARM_PARMHOLDER_H
#include <map>
#include "ParmEnum.h"
#include "../TypeNameHolder.h"
namespace Cpptraj {
namespace Parm {
/// Used to associate atom type names with an object (parameter etc)
template <class T> class ParmHolder {
    typedef std::pair<TypeNameHolder,T> Bpair;
    typedef std::map<TypeNameHolder,T> Bmap;
  public:
    /// CONSTRUCTOR
    ParmHolder() {}
    /// Clear all parameters
    void clear()              { bpmap_.clear(); }
    /// \return Number of parameters
    size_t size()       const { return bpmap_.size(); }
    /// \return true if no parameters
    bool empty()        const { return bpmap_.empty(); }
    /// \return Last parameter to be overwritten from AddParm()
    T const& PreviousParm() const { return previousParm_; }
    /// Add (or update if allowed) given parameter to holder.
    RetType AddParm(TypeNameHolder const& typesIn, T const& bp, bool allowUpdate) {
      // Ensure types are sorted
      TypeNameHolder types = typesIn;
      types.SortNames();
      // Check if parm for these types exist
      typename Bmap::iterator it = bpmap_.lower_bound( types );
      if (it == bpmap_.end() || it->first != types) {
        // New parm
        it = bpmap_.insert(it, Bpair(types, bp));
      } else {
        if (bp < it->second || it->second < bp) {
          // Potential update of existing parameter
          if (allowUpdate) {
            previousParm_ = it->second;
            it->second = bp;
            return UPDATED;
          } else {
            return ERR;
          }
        } else {
          //mprintf("DEBUG: Existing parameter:");
          //for (TypeNameHolder::const_iterator it = types.begin(); it != types.end(); ++it)
          //  mprintf(" '%s'", *(*it));
          //mprintf("\n");
          return SAME;
        }
      }
      return ADDED;
    }
    /// Constant iterator
    typedef typename Bmap::const_iterator const_iterator;
    /// \return constant iterator to beginning
    const_iterator begin() const { return bpmap_.begin(); }
    /// \return constant iterator to end.
    const_iterator end()   const { return bpmap_.end();   }
    /// Iterator
    typedef typename Bmap::iterator iterator;
    /// \return iterator to beginning
    iterator begin() { return bpmap_.begin(); }
    /// \return iterator to end
    iterator end()   { return bpmap_.end();   }
    /// \return Parameter matching given types, or empty parameter if not found.
    T FindParam(TypeNameHolder const& typesIn, bool& found) const {
      TypeNameHolder types = typesIn;
      types.SortNames();
      typename Bmap::const_iterator it = bpmap_.find( types );
      if (it == bpmap_.end()) {
        found = false;
        return T();
      }
      found = true;
      return it->second;
    }
    /// \return iterator to parameter matching the given types.
    iterator GetParam(TypeNameHolder const& typesIn) {
      TypeNameHolder types = typesIn;
      types.SortNames();
      return bpmap_.find( types );
    }
    /// \return const iterator to parameter matching the given types.
    const_iterator GetParam(TypeNameHolder const& typesIn) const {
      TypeNameHolder types = typesIn;
      types.SortNames();
      return bpmap_.find( types );
    }
    /// \return size in memory in bytes
    size_t DataSize() const {
      if (bpmap_.empty()) return 0;
      const_iterator elt0 = begin();
      // Assume all TypeNameHolders are the same size
      return (bpmap_.size() * elt0->first.DataSize()) +
             (bpmap_.size() * sizeof(T)) +
             sizeof(Bmap);
    }
  private:
    Bmap bpmap_;
    T previousParm_; ///< When parameter is updated, store previous value.
};
} // END namespace Parm
} // END namespace Cpptraj
#endif

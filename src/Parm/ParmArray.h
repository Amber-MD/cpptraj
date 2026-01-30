#ifndef INC_PARM_PARMARRAY_H
#define INC_PARM_PARMARRAY_H
#include <vector>
#include <utility> // std::pair
#include "ParmEnum.h"
#include "../TypeNameHolder.h"
//#incl ude "CpptrajStdio.h" // DEBUG
namespace Cpptraj {
namespace Parm {
/// Used to associate atom type names with an object (parameter etc) in an array
template <class T> class ParmArray {
    typedef std::pair<TypeNameHolder,T> Bpair;
    typedef std::vector<Bpair> Bmap;
  public:
    /// CONSTRUCTOR
    ParmArray() {}
    /// Clear all parameters
    void clear()              { bpmap_.clear(); }
    /// \return Number of parameters
    size_t size()       const { return bpmap_.size(); }
    /// \return true if no parameters
    bool empty()        const { return bpmap_.empty(); }
    /// \return Last parameter to be overwritten from AddParm()
    T const& PreviousParm() const { return previousParm_; }
    /// Add (or update if allowed) given parameter to holder.
    RetType AddParm(TypeNameHolder const& types, T const& bp, bool allowUpdate) {
      // Check if parm for these types exist
      typename Bmap::iterator it = bpmap_.begin();
      for (; it != bpmap_.end(); ++it)
        if (it->first.Match_NoWC( types )) break;
      if (it == bpmap_.end()) {
        // New parm
        //mprintf("DEBUG: New parameter:");
        //for (TypeNameHolder::const_iterator it = types.begin(); it != types.end(); ++it)
        //  mprintf(" '%s'", *(*it));
        //mprintf("\n");
        bpmap_.push_back( Bpair(types, bp) );
      } else {
        if (bp < it->second || it->second < bp) {
          //mprintf("DEBUG: Potential update of existing parameter:");
          //for (TypeNameHolder::const_iterator it = types.begin(); it != types.end(); ++it)
          //  mprintf(" '%s'", *(*it));
          //mprintf("\n");
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
    /// \return const iterator to parameter matching the given types.
    const_iterator GetParam(TypeNameHolder const& types) const {
      for (const_iterator it = bpmap_.begin(); it != bpmap_.end(); ++it)
        if (it->first.Match_NoWC( types )) return it;
      return bpmap_.end();
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

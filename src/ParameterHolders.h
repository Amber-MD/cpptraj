#ifndef INC_PARAMETERHOLDERS_H
#define INC_PARAMETERHOLDERS_H
#include <vector>
#include <utility> // std::pair
#include "TypeNameHolder.h"
#include "ParameterTypes.h"
//#incl ude "CpptrajStdio.h" // DEBUG

namespace ParameterHolders {
  enum RetType { ADDED = 0, SAME, UPDATED, ERR };
} /* END namespace ParameterHolders */

/// Used to associate atom type names with an object (parameter etc)
template <class T> class ParmHolder {
    // TODO may want to actually use a map one day for performance reasons.
    typedef std::pair<TypeNameHolder,T> Bpair;
    typedef std::vector<Bpair> Bmap;
  public:
    ParmHolder() {}
    void clear()              { bpmap_.clear(); }
    size_t size()       const { return bpmap_.size(); }
    bool empty()        const { return bpmap_.empty(); }
    /// Add (or update if allowed) given parameter to holder.
    ParameterHolders::RetType AddParm(TypeNameHolder const& types, T const& bp, bool allowUpdate) {
      // Check if parm for these types exist
      typename Bmap::iterator it = bpmap_.begin();
      for (; it != bpmap_.end(); ++it)
        if (it->first == types) break;
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
            it->second = bp;
            return ParameterHolders::UPDATED;
          } else {
            return ParameterHolders::ERR;
          }
        } else {
          //mprintf("DEBUG: Existing parameter:");
          //for (TypeNameHolder::const_iterator it = types.begin(); it != types.end(); ++it)
          //  mprintf(" '%s'", *(*it));
          //mprintf("\n");
          return ParameterHolders::SAME;
        }
      }
      return ParameterHolders::ADDED;
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
    T FindParam(TypeNameHolder const& types, bool& found) const { // TODO only use GetParam()?
      found = true;
      for (const_iterator it = begin(); it != end(); ++it)
        if (it->first == types) return it->second;
      found = false;
      return T();
    }
    /// \return iterator to parameter matching the given types.
    iterator GetParam(TypeNameHolder const& types) {
      for (iterator it = bpmap_.begin(); it != bpmap_.end(); ++it)
        if (it->first == types) return it;
      return bpmap_.end();
    }
    /// \return const iterator to parameter matching the given types.
    const_iterator GetParam(TypeNameHolder const& types) const {
      for (const_iterator it = bpmap_.begin(); it != bpmap_.end(); ++it)
        if (it->first == types) return it;
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
};

// -----------------------------------------------------------------------------
/// Specialized class for associating atom types with dihedral parameters.
/** NOTE: Instead of using a specialize template here I'm creating a new
  *       class because while I want AddParm() to accept DihedralParmType,
  *       I want FindParam to return an array of DihedralParmType, one for
  *       each unique multiplicity.
  */
class DihedralParmHolder {
    typedef std::pair<TypeNameHolder,DihedralParmArray> Bpair;
    typedef std::vector<Bpair> Bmap;
  public:
    DihedralParmHolder() {}
    void clear()              { bpmap_.clear();        }
    size_t size()       const { return bpmap_.size();  }
    bool empty()        const { return bpmap_.empty(); }
    /** Add (or update) a single dihedral parameter for given atom types. */
    ParameterHolders::RetType
    AddParm(TypeNameHolder const& types, DihedralParmType const& dp, bool allowUpdate) {
      // Check if parm for these types exist
      Bmap::iterator it0 = bpmap_.begin();
      for (; it0 != bpmap_.end(); ++it0)
      {
        if (it0->first == types)
          break;
      }
      if (it0 == bpmap_.end()) {
        // Brand new dihedral for these types.
        //mprintf("DEBUG: New dihedral parm: %s %s %s %s pk=%12.4f pn=%12.4f pp=%12.4f\n",
        //        *types[0], *types[1], *types[2], *types[3], dp.Pk(), dp.Pn(), dp.Phase());
        bpmap_.push_back( Bpair(types, DihedralParmArray(1, dp)) );
      } else {
        // If we are here types match - check multiplicity.
        DihedralParmArray::iterator it1 = it0->second.begin();
        for (; it1 != it0->second.end(); ++it1)
        {
          if (it1->Pn() == dp.Pn())
            break;
        }
        if (it1 == it0->second.end()) {
          // Brand new multiplicity for this dihedral.
          //mprintf("DEBUG: Dihedral new mult: %s %s %s %s pk=%12.4f pn=%12.4f pp=%12.4f\n",
          //        *types[0], *types[1], *types[2], *types[3], dp.Pk(), dp.Pn(), dp.Phase());
          it0->second.push_back( dp );
        } else {
          if (dp < *it1 || *it1 < dp) {
            //mprintf("DEBUG: Attempt dihedral update mult (allow=%i): %s %s %s %s pk=%6.2f pn=%3.1f pp=%6.3f (orig pk=%6.2f pn=%3.1f pp=%6.3f )\n",
            //        (int)allowUpdate, *types[0], *types[1], *types[2], *types[3], dp.Pk(), dp.Pn(), dp.Phase(), it1->Pk(), it1->Pn(), it1->Phase());
            if (allowUpdate) {
              *it1 = dp;
              return ParameterHolders::UPDATED;
            } else {
              return ParameterHolders::ERR;
            }
          } else
            return ParameterHolders::SAME;
        }
      }
      return ParameterHolders::ADDED;
    }

    /** This version takes an array of dihedral parameters. */
    ParameterHolders::RetType
    AddParm(TypeNameHolder const& types, DihedralParmArray const& dpa, bool allowUpdate) {
      // Check if parm for these types exist
      Bmap::iterator it0 = bpmap_.begin();
      for (; it0 != bpmap_.end(); ++it0)
      {
        if (it0->first == types)
          break;
      }
      if (it0 == bpmap_.end()) {
        // Brand new dihedral for these types.
        bpmap_.push_back( Bpair(types, dpa) );
      } else {
        if (!allowUpdate) return ParameterHolders::ERR;
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
          it0->second = dpa;
          return ParameterHolders::UPDATED;
        } else
          return ParameterHolders::SAME;
      }
      return ParameterHolders::ADDED;
    }

    typedef typename Bmap::const_iterator const_iterator;
    const_iterator begin() const { return bpmap_.begin(); }
    const_iterator end()   const { return bpmap_.end();   }
    /// \return Array of dihedral parameters matching given atom types.
    DihedralParmArray FindParam(TypeNameHolder const& types, bool& found) const {
      found = true;
      for (const_iterator it = begin(); it != end(); ++it)
        if (it->first == types) return it->second;
      found = false;
      return DihedralParmArray();
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
  private:
    Bmap bpmap_;
};
#endif

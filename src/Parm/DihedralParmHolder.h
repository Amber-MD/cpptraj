#ifndef INC_PARM_DIHEDRALPARMHOLDER_H
#define INC_PARM_DIHEDRALPARMHOLDER_H
#include <vector>
#include <utility> // std::pair
#include "ParmEnum.h"
#include "../ParameterTypes.h"
#include "../TypeNameHolder.h"
//#incl ude "CpptrajStdio.h" // DEBUG
namespace Cpptraj {
namespace Parm {
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
    DihedralParmType const& PreviousParm() const { return previousParm_; }
    DihedralParmArray const& PreviousArray() const { return previousArray_; }
    /// Set wildcard character
    void SetWildcard(char wc) { wc_ = NameType(std::string(1, wc)); }
    /// Add (or update) a single dihedral parameter for given atom types.
    RetType AddParm(TypeNameHolder const& types, DihedralParmType const& dp, bool allowUpdate) {
      // Check if parm for these types exist
      Bmap::iterator it0 = bpmap_.begin();
      for (; it0 != bpmap_.end(); ++it0)
      {
        if (it0->first.Match_NoWC( types ))
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
          if (FEQ(it1->Pn(), dp.Pn()))
            break;
        }
        if (it1 == it0->second.end()) {
          // Brand new multiplicity for this dihedral.
          //mprintf("DEBUG: Dihedral new mult: %s %s %s %s pk=%12.4f pn=%12.4f pp=%12.4f\n",
          //        *types[0], *types[1], *types[2], *types[3], dp.Pk(), dp.Pn(), dp.Phase());
          if (it0->second.empty())
            it0->second.push_back( dp );
          else if (dp.Pn() > it0->second.back().Pn())
            it0->second.push_back( dp );
          else {
            // Try to keep multiplicities in order.
            DihedralParmArray sorted;
            bool isInserted = false;
            for (DihedralParmArray::const_iterator jt = it0->second.begin(); jt != it0->second.end(); ++jt) {
              if (!isInserted) {
                if (dp.Pn() < jt->Pn()) {
                  sorted.push_back( dp );
                  isInserted = true;
                }
              }
              sorted.push_back( *jt );
            }
            it0->second = sorted;
          }
        } else {
          if (dp < *it1 || *it1 < dp) {
            //mprintf("DEBUG: Attempt dihedral update mult (allow=%i): %s %s %s %s pk=%6.2f pn=%3.1f pp=%6.3f (orig pk=%6.2f pn=%3.1f pp=%6.3f )\n",
            //        (int)allowUpdate, *types[0], *types[1], *types[2], *types[3], dp.Pk(), dp.Pn(), dp.Phase(), it1->Pk(), it1->Pn(), it1->Phase());
            if (allowUpdate) {
              previousParm_ = *it1;
              *it1 = dp;
              return UPDATED;
            } else {
              return ERR;
            }
          } else
            return SAME;
        }
      }
      return ADDED;
    }
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
          previousArray_ = it0->second;
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
    DihedralParmType previousParm_; ///< When parameter is updated, store previous value. TODO remove?
    DihedralParmArray previousArray_; ///< When array is updated, store previous value
};
} // END namespace Parm
} // END namespace Cpptraj
#endif

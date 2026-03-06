#ifndef INC_PARM_IMPROPERPARMHOLDER_H
#define INC_PARM_IMPROPERPARMHOLDER_H
#include "DihedralParmHolder.h"
//#incl ude "CpptrajStdio.h" // DEBUG
namespace Cpptraj {
namespace Parm {
// -----------------------------------------------------------------------------
/// Specialized class for associating atom types with improper parameters.
/** Impropers are a little tricky in that by convention the third atom is
  * the central atom, and all other atoms can be in any order.
  * The Amber convention is usually (but not always) to have the non-central
  * improper atom types sorted alphabetically, with wildcards given
  * precedence, but this is not always the case and does not always work.
  * For example, using straight up backwards/forwards matching, the wildcard
  * type X-X-CW-H4 will not match the given alphabetized type C*-H4-CW-NA.
  * All combinations of A1, A2, and A4 should be checked.
  */
class ImproperParmHolder : private DihedralParmHolder {
    /// Function for matching wildcards (WildCard Match)
    static inline bool wcm(NameType const& t0, NameType const& t1, NameType const& wc) {
      return (t0 == wc || t0 == t1);
    }
    /// Function for swapping integers
    static inline void swp(int& i1, int& i2) { int tmp = i1; i1 = i2; i2 = tmp; }
  public:
    /// Denote returned parameter atom type order
    enum OrderType { O_013,
                     O_031,
                     O_103,
                     O_130,
                     O_301,
                     O_310 };
    ImproperParmHolder() : require_exact_match_(false) {}
    /// \return Number of improper parameter sets
    size_t size()              const { return DihedralParmHolder::size();  }
    /// \return True if no parameters
    bool empty()               const { return DihedralParmHolder::empty(); }
    /// \return Last parameter to be overwritten from AddParm()
    DihedralParmArray const& PreviousParm() const { return DihedralParmHolder::PreviousParm(); }
    /// \return Wildcard
    NameType const& Wildcard() const { return wc_; }
    /// \return True if an exact match is required to find a parameter
    bool RequireExactMatch()   const { return require_exact_match_; }
    /// const iterator
    typedef typename DihedralParmHolder::const_iterator const_iterator;
    /// const iterator to beginning of parameters
    const_iterator begin()     const { return DihedralParmHolder::begin(); }
    /// const iterator to end of parameters
    const_iterator end()       const { return DihedralParmHolder::end();   }
    /// Set Wildcard char
    void SetWildcard(char wc) { DihedralParmHolder::SetWildcard(wc); }
    /// Set Wildcard
    void SetWildcard(NameType const& wc) { wc_ = wc; }
    /// Indicate whether exact type matches are required to find parameters
    void SetRequireExactMatch(bool b) { require_exact_match_ = b; }
/*
    /// Add (or update) a single improper parameter for given atom types.
    RetType AddParm(TypeNameHolder const& types, DihedralParmType const& dp, bool allowUpdate) {
      return DihedralParmHolder::AddParm( types, dp, allowUpdate );
    }
*/
    /// Add array of improper parameters, one for each multiplicity.
    RetType AddParm(TypeNameHolder const& types, DihedralParmArray const& dpa, bool allowUpdate) {
      return DihedralParmHolder::AddParm( types, dpa, allowUpdate );
    }
  private:
    /// Remap the given improper according to the desired order, correct for wildcards in given types
    void ReorderImproper(DihedralType& imp, OrderType order, TypeNameHolder const& types) const {
      switch (order) {
        case O_013 : break;
        case O_031 : swp( imp.ChangeA2(), imp.ChangeA4() ); break;
        case O_103 : swp( imp.ChangeA1(), imp.ChangeA2() ); break;
        case O_130 : swp( imp.ChangeA1(), imp.ChangeA4() ); swp( imp.ChangeA1(), imp.ChangeA2() ); break;
        case O_301 : swp( imp.ChangeA2(), imp.ChangeA4() ); swp( imp.ChangeA1(), imp.ChangeA2() ); break;
        case O_310 : swp( imp.ChangeA1(), imp.ChangeA4() ); break;
      }
      if (wc_.len() > 0 && types.Size() > 0) {
        // If there are wildcards, need to order matching types by atom index.
        // General order of imp should now match types due to above swaps.
        bool wc1 = (types[0] == wc_);
        bool wc2 = (types[1] == wc_);
        bool wc4 = (types[3] == wc_);
        //mprintf("DEBUG: %s WC0=%i %s WC1=%i %s WC3=%i\n", *types[0], (int)wc1, *types[1], (int)wc2, *types[3], (int)wc4);
        if (wc4 && wc2 && wc1) {
          // All three wildcard - should be rare but need to check.
          if (imp.A1() > imp.A2()) swp(imp.ChangeA1(), imp.ChangeA2());
          if (imp.A2() > imp.A4()) swp(imp.ChangeA2(), imp.ChangeA4());
          if (imp.A1() > imp.A2()) swp(imp.ChangeA1(), imp.ChangeA2());
          if (imp.A2() > imp.A4()) swp(imp.ChangeA2(), imp.ChangeA4());
        } else if (wc1 && wc2) {
          if (imp.A1() > imp.A2()) swp(imp.ChangeA1(), imp.ChangeA2());
        } else if (wc1 && wc4) {
          if (imp.A1() > imp.A4()) swp(imp.ChangeA1(), imp.ChangeA4());
        } else if (wc2 && wc4) {
          if (imp.A2() > imp.A4()) swp(imp.ChangeA2(), imp.ChangeA4());
        }
      }
    }
  public:
    /// Remap the given improper according to the desired order. Used by unit test to check that reordering works.
    void ReorderImproper(DihedralType& imp, OrderType order) const {
      ReorderImproper(imp, order, TypeNameHolder());
    }
    /// Get improper parameters matching given atom types. If found, improper will be reordered to match parameter order.
    const_iterator GetParam(TypeNameHolder const& types, DihedralType& imp, bool& reordered) const {
      //mprintf("DEBUG: FindParam wc=%s Inco=%s-%s-%s-%s\n",*wc_, *(types[0]), *(types[1]),   *(types[2]),   *(types[3]));
      reordered = false;
      OrderType lastOrder_ = O_013;
      const_iterator it = begin();
      // If we require an exact match, look for that first.
      if (require_exact_match_) {
        for (; it != end(); ++it) {
          TypeNameHolder const& myTypes = it->first;
          if (myTypes[2] == types[2] && myTypes[0] == types[0] &&
              myTypes[1] == types[1] && myTypes[3] == types[3])
          {
            break;
          }
        }
      } else {
        // Inexact match. First, no wildcard
        for (; it != end(); ++it) {
          TypeNameHolder const& myTypes = it->first;
          // Central (third) type must match
          if (myTypes[2] == types[2]) {
            //mprintf("DEBUG: FindParam (improper) central atom match %s", *(types[2]));
            //mprintf(" This=%s-%s-%s-%s", *(myTypes[0]), *(myTypes[1]), *(myTypes[2]), *(myTypes[3]));
            //mprintf(" Inco=%s-%s-%s-%s\n", *(types[0]), *(types[1]),   *(types[2]),   *(types[3]));
            // Try all permutations
            if (       myTypes[0] == types[0] && myTypes[1] == types[1] && myTypes[3] == types[3]) { // 0 1 2 3
                lastOrder_ = O_013; break;
            } else if (myTypes[0] == types[0] && myTypes[1] == types[3] && myTypes[3] == types[1]) { // 0 3 2 1
                lastOrder_ = O_031; break;
            } else if (myTypes[0] == types[1] && myTypes[1] == types[0] && myTypes[3] == types[3]) { // 1 0 2 3
                lastOrder_ = O_103; break;
            } else if (myTypes[0] == types[1] && myTypes[1] == types[3] && myTypes[3] == types[0]) { // 1 3 2 0
                lastOrder_ = O_130; break;
            } else if (myTypes[0] == types[3] && myTypes[1] == types[0] && myTypes[3] == types[1]) { // 3 0 2 1
                lastOrder_ = O_301; break;
            } else if (myTypes[0] == types[3] && myTypes[1] == types[1] && myTypes[3] == types[0]) { // 3 1 2 0
                lastOrder_ = O_310; break;
            }
          }
        } // END loop over parameters
        // Second, do wildcard matches if wildcard is set.
        if (it == end() && wc_.len() > 0) {
          it = begin();
          for (; it != end(); ++it) {
            TypeNameHolder const& myTypes = it->first;
            // Central (third) type must match
            if (wcm(myTypes[2], types[2], wc_)) {
              // Try all permutations
              if (       wcm(myTypes[0], types[0], wc_) && wcm(myTypes[1], types[1], wc_) && wcm(myTypes[3], types[3], wc_)) { // 0 1 2 3
                  lastOrder_ = O_013; break;
              } else if (wcm(myTypes[0], types[0], wc_) && wcm(myTypes[1], types[3], wc_) && wcm(myTypes[3], types[1], wc_)) { // 0 3 2 1
                  lastOrder_ = O_031; break;
              } else if (wcm(myTypes[0], types[1], wc_) && wcm(myTypes[1], types[0], wc_) && wcm(myTypes[3], types[3], wc_)) { // 1 0 2 3
                  lastOrder_ = O_103; break;
              } else if (wcm(myTypes[0], types[1], wc_) && wcm(myTypes[1], types[3], wc_) && wcm(myTypes[3], types[0], wc_)) { // 1 3 2 0
                  lastOrder_ = O_130; break;
              } else if (wcm(myTypes[0], types[3], wc_) && wcm(myTypes[1], types[0], wc_) && wcm(myTypes[3], types[1], wc_)) { // 3 0 2 1
                  lastOrder_ = O_301; break;
              } else if (wcm(myTypes[0], types[3], wc_) && wcm(myTypes[1], types[1], wc_) && wcm(myTypes[3], types[0], wc_)) { // 3 1 2 0
                  lastOrder_ = O_310; break;
              }
            }
          } // END loop over parameters
        } // END wildcard matches
      } // END require exact match
      if (it != end()) {
        // We have found a parameter. Do any reordering.
        if (lastOrder_ != O_013) reordered = true;
        ReorderImproper( imp, lastOrder_, it->first );
      }
      return it;
    } // END GetParam()
    /// Get improper parameters matching given atom types. If found, improper will be reordered to match parameter order.
    const_iterator GetParam(TypeNameHolder const& types) const {
      DihedralType imp;
      bool reordered = false;
      return GetParam( types, imp, reordered);
    }
    /// \return Array of improper parameters matching given atom types. Improper will be reordered to match parameter order.
    DihedralParmArray FindParam(TypeNameHolder const& types, bool& found, DihedralType& imp, bool& reordered) const {
      const_iterator it = GetParam( types, imp, reordered );
      if (it == end()) {
        found = false;
        return DihedralParmArray();
      } else {
        found = true;
        return it->second;
      }
    } // END FindParam()
    /// \return Dihedral parm array corresponding to types. Use by unit test.
    DihedralParmArray FindParam(TypeNameHolder const& types, bool& found) const {
      DihedralType blank;
      bool reordered;
      return FindParam(types, found, blank, reordered);
    }
    /// \return size in memory in bytes
    size_t DataSize() const { return DihedralParmHolder::DataSize(); }
  private:
    bool require_exact_match_; ///< If true, types must match in exact order when finding parameters
};
} // END namespace Parm
} // END namespace Cpptraj
#endif

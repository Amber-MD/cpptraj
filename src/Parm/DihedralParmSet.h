#ifndef INC_PARM_DIHEDRALPARMSET_H
#define INC_PARM_DIHEDRALPARMSET_H
//#include <map>
#include "ParmEnum.h"
#include "../ParameterTypes.h"
#include "../TypeNameHolder.h"
namespace Cpptraj {
namespace Parm {
class DihedralParmHolder;
class ImproperParmHolder;
/// Used for reading in dihedral parameters.
/** When dihedral parameters are read in one at a time it is not
  * guaranteed that parameters for all multiplcities will be
  * present. This allows dihedral parameters to be read in
  * one multiplicit at a time and then be checked for
  * consistency afterwards.
  */
class DihedralParmSet {
    typedef std::pair<TypeNameHolder, DihedralParmArray> Dpair;
    //typedef std::map<TypeNameHolder, DihedralParmArray> Dmap;
    typedef std::vector<Dpair> Dmap;
  public:
    /// CONSTRUCTOR
    DihedralParmSet() : debug_(0) {}
    /// CONSTRUCTOR - with debug level
    DihedralParmSet(int debugIn) : debug_(debugIn) {}
    /// \return Last parameter to be overwritten from AddParm()
    DihedralParmType const& PreviousParm() const { return previousParm_; }

    typedef Dmap::const_iterator const_iterator;
    const_iterator begin() const { return dihparm_.begin(); }
    const_iterator end()   const { return dihparm_.end(); }

    /// Put all dihedral parameters into given dihedral parameter holder.
    int ToDihParm(DihedralParmHolder&) const;
    /// Put all improper parameters into given improper parameter holder.
    int ToImpParm(ImproperParmHolder&) const;

    /// Add (or update) a single dihedral or improper parameter for given atom types.
    RetType AddDihParm(TypeNameHolder const& types, DihedralParmType const& dp, bool allowUpdate) {
      // Ensure types are sorted
      //TypeNameHolder types = typesIn;
      //types.SortNames();
      // Check if dihedral parm for these types exist
      //Dmap::iterator it0 = dihparm_.lower_bound( types );
      //if (it0 == dihparm_.end() || it0->first != types) {
      //  // Brand new dihedral for these types
      //  it0 = dihparm_.insert(it0, Dpair(types, DihedralParmArray(1, dp)));
      // Check if parm for these types exist
      Dmap::iterator it0 = dihparm_.begin();
      for (; it0 != dihparm_.end(); ++it0)
      {
        if (it0->first.Match_NoWC( types ))
          break;
      }
      if (it0 == dihparm_.end()) {
        // Brand new dihedral for these types.
        //mprintf("DEBUG: New dihedral parm: %s %s %s %s pk=%12.4f pn=%12.4f pp=%12.4f\n",
        //        *types[0], *types[1], *types[2], *types[3], dp.Pk(), dp.Pn(), dp.Phase());
        dihparm_.push_back( Dpair(types, DihedralParmArray(1, dp)) );
      } else {
        // If we are here types match - check multiplicity
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
          // Multiplicity already exists for this dihedral.
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
  private:
    static inline void check_mult(std::string const&, TypeNameHolder const&, DihedralParmArray const&);

    Dmap dihparm_;
    DihedralParmType previousParm_; ///< When parameter is updated, store previous value.
    int debug_; ///< Debug level
};
}
}
#endif

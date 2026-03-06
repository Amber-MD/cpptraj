#ifndef INC_PARM_MERGE_H
#define INC_PARM_MERGE_H
#include <vector>
class Atom;
class Residue;
class BondArray;
class BondParmArray;
class AngleArray;
class AngleParmArray;
class DihedralArray;
class DihedralParmArray;
class CmapArray;
class CmapGridArray;
namespace Cpptraj {
namespace Parm {
/// Routines used to merge Topology files.
/** Compile with -DCPPTRAJ_DEBUG_MERGE for extra debug info. */
    /// Shorthand for array of Atoms
    typedef std::vector<Atom> AtArray;
    /// Shorthand for array of residues
    typedef std::vector<Residue> ResArray;
    /// Merge two pairs of bond arrays and their parameters.
    void MergeBondArrays(bool, BondArray&, BondArray&, BondParmArray&, AtArray const&,
                         BondArray const&, BondArray const&, BondParmArray const&, AtArray const&);
    /// Merge two pairs of angle arrays and their parameters.
    void MergeAngleArrays(bool, AngleArray&, AngleArray&, AngleParmArray&, AtArray const&,
                          AngleArray const&, AngleArray const&, AngleParmArray const&, AtArray const&);
    /// Merge two pairs of dihedral arrays and their parameters.
    void MergeDihedralArrays(DihedralArray&, DihedralArray&, DihedralParmArray&, AtArray const&,
                             DihedralArray const&, DihedralArray const&, DihedralParmArray const&, AtArray const&);
    /// Merge cmap arrays and their parameters
    void MergeCmapArrays(CmapArray&, CmapGridArray&, AtArray const&, ResArray const&,
                         CmapArray const&, CmapGridArray const&, AtArray const&, ResArray const&);
    /// Merge two bond arrays and their parameters
    void MergeBondArray(BondArray&, BondParmArray&, AtArray const&,
                        BondArray const&, BondParmArray const&, AtArray const&);
    /// Merge two improper arrays and their parameters
    void MergeImproperArray(DihedralArray&, DihedralParmArray&, AtArray const&,
                            DihedralArray const&, DihedralParmArray const&, AtArray const&);

// -----------------------------------------------------------------------------
/** This template can be used when doing Append() on a generic std::vector array
  * of type T. The array will be appended to a given array of the same type.
  * If one is empty and the other is not, values will be filled in if necessary.
  */
template <class T> class TopVecAppend {
  public:
    /// CONSTRUCTOR
    TopVecAppend() {}
    /// Append current array to given array of same type
    void Append(std::vector<T>& arrayOut, std::vector<T> const& arrayToAdd, unsigned int expectedSize)
    {
      if (arrayToAdd.empty() && arrayOut.empty()) {
        // Both arrays are empty. Nothing to do.
        return;
      } else if (arrayToAdd.empty()) {
        // The current array is empty but the given array is not. Fill in 
        // array to append with blank values.
        for (unsigned int idx = 0; idx != expectedSize; idx++)
          arrayOut.push_back( T() );
      } else {
        // Append current array to array to given array. TODO use std::copy?
        for (typename std::vector<T>::const_iterator it = arrayToAdd.begin(); it != arrayToAdd.end(); ++it)
          arrayOut.push_back( *it );
      }
    }
};
}
}
#endif

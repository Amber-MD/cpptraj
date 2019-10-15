#ifndef INC_DATASET_TENSOR_H
#define INC_DATASET_TENSOR_H
#include <vector>
#include "SymmetricTensor.h"
#include "DataSet.h"
/// Hold an array of 3x3 symmetric tensor values. The array can be sparse.
/** NOTE: We use a floating point type array for indices in case we want to
  * hold time values at some point.
  */
class DataSet_Tensor : public DataSet {
    typedef SymmetricTensor<double> Ttype;
    typedef std::vector<Ttype> Tarray;
    typedef std::vector<double> Farray;
  public:
    /// CONSTRUCTOR
    DataSet_Tensor();
    /// Allocator
    static DataSet* Alloc() { return (DataSet*)new DataSet_Tensor(); }

    // ----- DataSet functions -------------------
    size_t Size() const { return Data_.size(); }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&);
#   endif
    void Info()                    const { return;                    }
    int Allocate(SizeArray const&);
    void Add( size_t, const void* );
    void WriteBuffer(CpptrajFile&, SizeArray const&) const;
    double Coord(unsigned int, size_t) const;
    int Append(DataSet*);
    size_t MemUsageInBytes() const;
  private:
    Farray Xvals_; /// < Hold the index values
    Tarray Data_;  ///< Hold the tensor values
};
#endif

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
  public:
    typedef SymmetricTensor<double> Ttype;

    /// CONSTRUCTOR
    DataSet_Tensor();
    /// Allocator
    static DataSet* Alloc() { return (DataSet*)new DataSet_Tensor(); }

    /// \return Index at specified position
    double Xvals(unsigned int idx) const { return Xvals_[idx]; }
    /// \return Tensor at specified position
    Ttype Tensor(unsigned int idx) const { return Data_[idx]; }

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
    typedef std::vector<Ttype> Tarray;
    typedef std::vector<double> Farray;

    Farray Xvals_; /// < Hold the index values
    Tarray Data_;  ///< Hold the tensor values
};
#endif

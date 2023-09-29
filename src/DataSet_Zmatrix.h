#ifndef INC_DATASET_ZMATRIX_H
#define INC_DATASET_ZMATRIX_H
#include "DataSet.h"
namespace Cpptraj {
namespace Structure {
class Zmatrix;
}
}
/// Hold a Z matrix 
class DataSet_Zmatrix : public DataSet {
  public:
    DataSet_Zmatrix();
    ~DataSet_Zmatrix();
    static DataSet* Alloc() { return (DataSet*)new DataSet_Zmatrix(); }
    // ----- DataSet functions -------------------
    size_t Size() const;
    int Allocate(SizeArray const&);
    void WriteBuffer(CpptrajFile&, SizeArray const&) const;
    size_t MemUsageInBytes() const;
    void Info()                                      const { return; }
    void Add(size_t, const void*)                          { return; }
    int Append(DataSet*)                                   { return 1; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    // -------------------------------------------
    Cpptraj::Structure::Zmatrix* Zptr() const { return zmatrix_; }
  private:
    Cpptraj::Structure::Zmatrix* zmatrix_;
};
#endif

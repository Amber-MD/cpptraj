#ifndef INC_DATASET_PARAMETERS_H
#define INC_DATASET_PARAMETERS_H
#include "DataSet.h"
/// <Enter description of DataSet_Parameters here>
class DataSet_Parameters : public DataSet {
  public:
    DataSet_Parameters();
    static DataSet* Alloc() { return (DataSet*)new DataSet_Parameters(); }
    // ----- DataSet functions -------------------
    size_t Size()                                    const { return 0; }
    void Info()                                      const { return; }
    int Allocate(SizeArray const&)                         { return 1; }
    void Add(size_t, const void*)                          { return; }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; }
    int Append(DataSet*)                                   { return 1; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    // -------------------------------------------
  private:
    
};
#endif

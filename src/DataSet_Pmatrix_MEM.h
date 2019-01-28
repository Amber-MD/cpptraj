#ifndef INC_DATASET_PMATRIX_MEM_H
#define INC_DATASET_PMATRIX_MEM_H
#include "DataSet.h"
#include "Cluster/PairwiseMatrix_MEM.h"
/// Cluster pairwise matrix, distances stored in memory.
class DataSet_Pmatrix_MEM : public DataSet, Cpptraj::Cluster::PairwiseMatrix_MEM {
  public:
    DataSet_Pmatrix_MEM() : DataSet(PMATRIX_MEM, CLUSTERMATRIX,
                                    TextFormat(TextFormat::DOUBLE, 12, 4), 2) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_Pmatrix_MEM(); }
    // ----- DataSet functions -------------------
    size_t Size() const { return 0; } // TODO
    void Info()   const { return; }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; } // TODO
    int Allocate(SizeArray const&) { return 1; } // TODO
    void Add(size_t, const void*) {}
    int Append(DataSet*) { return 1; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    // -------------------------------------------
};
#endif

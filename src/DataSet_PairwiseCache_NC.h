#ifndef INC_DATASET_PAIRWISECACHE_NC_H
#define INC_DATASET_PAIRWISECACHE_NC_H
#include "DataSet_PairwiseCache.h"
#include "Cluster/Cmatrix_NC.h"
/// Cache pairwise distances in NetCDF file. 
class DataSet_PairwiseCache_NC : public DataSet_PairwiseCache {
  public:
    DataSet_PairwiseCache_NC() : DataSet_PairwiseCache(PMATRIX_NC) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_PairwiseCache_NC(); }
    // ----- DataSet functions -------------------
    size_t Size()                                    const { return file_.MatrixSize(); }
    void Info()                                      const { return; }
    int Allocate(SizeArray const&)                         { return 0; }
    void Add(size_t, const void*)                          { return; }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; }
    size_t MemUsageInBytes()                         const { return 0; }
//    int Append(DataSet*)                                   { return 1; }
//#   ifdef MPI
//    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
//#   endif
    // ----- PairwiseCache functions ------------- 
    double CachedDistance(unsigned int i1, unsigned int i2) const {
      return file_.GetCmatrixElement(i1, i2);
    }
    size_t Nrows() const { return file_.MatrixRows(); }
    void SetElement(int x, int y, double val) { file_.WriteCmatrixElement(x, y, val); }
    int SetupCache(unsigned int, Cframes const&, int, std::string const&);
    // -------------------------------------------
  private:
    Cpptraj::Cluster::Cmatrix_NC file_;  ///< Hold cached distances on disk
};
#endif

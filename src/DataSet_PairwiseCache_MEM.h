#ifndef INC_DATASET_PAIRWISECACHE_MEM_H
#define INC_DATASET_PAIRWISECACHE_MEM_H
#include "DataSet_PairwiseCache.h"
#include "Matrix.h"
/// <Enter description of DataSet_PairwiseCache_MEM here>
class DataSet_PairwiseCache_MEM : public DataSet_PairwiseCache {
  public:
    DataSet_PairwiseCache_MEM() : DataSet_PairwiseCache(PMATRIX_MEM) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_PairwiseCache_MEM(); }
    // ----- DataSet functions -------------------
    size_t Size()                                    const { return Mat_.size(); }
    void Info()                                      const { return; }
    int Allocate(SizeArray const&);
    void Add(size_t, const void*)                          { return; }
    void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; }
    size_t MemUsageInBytes()                         const { return Mat_.DataSize(); }
//    int Append(DataSet*)                                   { return 1; }
//#   ifdef MPI
//    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
//#   endif
    // -------------------------------------------
    //double GetFdist(int f1, int f2) const { return Mat_.element(FrameToMat()[f1], FrameToMat()[f2]); }
    //double Frame_Distance(int, int) const;
    // ----- PairwiseCache functions -------------
    double CachedDistance(unsigned int i1, unsigned int i2) const { return Mat_.element(i1, i2); }
    //int CacheDistances(Cframes const&);
    size_t Nrows() const { return Mat_.Nrows(); }
    void SetElement(int col, int row, double val) { Mat_.setElement(col, row, val); }
    int SetupCache(unsigned int, Cframes const&, int, std::string const&);
    // -------------------------------------------
    float*       Ptr()       { return Mat_.Ptr(); }
    float const* Ptr() const { return Mat_.Ptr(); }
  private:
    Matrix<float> Mat_;  ///< Hold cached distances
};
#endif

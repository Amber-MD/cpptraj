#ifndef INC_DATASET_PAIRWISECACHE_H
#define INC_DATASET_PAIRWISECACHE_H
#include "DataSet.h"
#include "Cluster/Cframes.h"
/// Used to cache distances in a potentially sparse square matrix 
class DataSet_PairwiseCache : public DataSet {
  public:
    DataSet_PairwiseCache() {}
    DataSet_PairwiseCache(DataType t) : DataSet(t, PWCACHE,
                                                TextFormat(TextFormat::DOUBLE, 12, 4), 2) {}
    virtual ~DataSet_PairwiseCache() {}

    typedef Cpptraj::Cluster::Cframes Cframes;
    //static DataSet* Alloc() { return (DataSet*)new DataSet_PairwiseCache(); }
    // ----- DataSet functions -------------------
    //size_t Size()                                    const { return 0; }
    //void Info()                                      const { return; }
    //int Allocate(SizeArray const&)                         { return 1; }
    //void Add(size_t, const void*)                          { return; }
    //void WriteBuffer(CpptrajFile&, SizeArray const&) const { return; }
    int Append(DataSet*)                                   { return 1; }
#   ifdef MPI
    int Sync(size_t, std::vector<int> const&, Parallel::Comm const&) { return 1; }
#   endif
    // -------------------------------------------
    /// \return distance between given cached frames (absolute frame numbers).
    //virtual double GetFdist(int, int) const = 0;
    /// \return distance between cached frames (internal indices).
    virtual double CachedDistance(unsigned int, unsigned int) const = 0;
    /// Print only cached distances.
    virtual void PrintCached() const = 0;
    /// Used to cache distances; expect internal indices, not absolute cluster frames.
    virtual void SetElement(int, int, double) = 0;
    /// Used to set up and allocate cache for total, frames to cache, sieve, and metric descrip.
    virtual int SetupCache(unsigned int, Cframes const&, int, std::string const&) = 0;
    // -------------------------------------------

    /// \return The array for converting frame numbers to internal cache indices.
    Cframes const& FrameToIdx()            const { return frameToIdx_; }
    /// \return True if given frames array matches cached frames.
    bool CachedFramesMatch(Cframes const&) const;
  protected:
    /// Set up frame to internal index array using given frames and total frames.
    int SetupFrameToIdx(Cframes const&, unsigned int);
  private:
    Cframes frameToIdx_; ///< Hold indices for all cached frames, -1 for non-cached.
};
#endif

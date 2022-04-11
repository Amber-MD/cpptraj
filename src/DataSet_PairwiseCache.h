#ifndef INC_DATASET_PAIRWISECACHE_H
#define INC_DATASET_PAIRWISECACHE_H
#include "DataSet.h"
#include "Cluster/Cframes.h"
/// Used to cache distances in a potentially sparse square symmetric matrix.
/** This is an abstract base class for storing pairwise distances
  * for cluster analysis. The frameToIdx_ array maps the matrix indices to
  * internal indices; indices for which there is no data are assigned -1.
  */
class DataSet_PairwiseCache : public DataSet {
  public:
    DataSet_PairwiseCache();
    /// This constructor should be used by child classes.
    DataSet_PairwiseCache(DataType);
    virtual ~DataSet_PairwiseCache() {}

    /// Class used to hold frame indices
    typedef Cpptraj::Cluster::Cframes Cframes;
    /// Class used to hold frame statuses (T = present)
    typedef std::vector<char> StatusArray;

    /// Character that means frame "present"; 'F' meaning frame was not sieve.
    static const char PRESENT_;
    /// Character that means frame "absent" - 'T' meaning frame was sieved.
    static const char ABSENT_; 

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
    /// \return Number of rows (i.e. number of frames cached) in pairwise cache. TODO DataSet_2D
    virtual size_t Nrows() const = 0;
    /// Used to cache distances; expect internal indices, not absolute cluster frames.
    virtual void SetElement(int, int, double) = 0;
    /// Used to set up and allocate cache for total, frames to cache, sieve, and metric descrip.
    virtual int SetupCache(unsigned int, Cframes const&, int, std::string const&) = 0;
    // -------------------------------------------

    /// \return The array for converting frame numbers to internal cache indices.
    Cframes const& FrameToIdx()            const { return frameToIdx_; }
    /// \return True if given frames array matches cached frames.
    bool CachedFramesMatch(Cframes const&) const;
    /// \return Internal sieve value for bookkeeping
    int SieveVal()                         const { return sieve_; }
    /// \return Internal metric description for bookkeeping
    std::string const& MetricDescrip()     const { return metricDescription_; }
    /// \return Array of frames present in the cache.
    Cframes PresentFrames()                const;
    /// Print only cached distances.
    void PrintCached() const;

    /// Used to set up frame indices from given status array. Record sieve value.
    int SetupFromStatus(StatusArray const&, int);
  protected:
    /// Set up frame to internal index array using given frames and total frames.
    int SetupFrameToIdx(Cframes const&, unsigned int);
    /// Set internal sieve value for bookkeeping (should call from SetupCache).
    void SetSieveVal(int s)                     { sieve_ = s;             }
    /// Set internal metric description for bookkeeping (should call from SetupCache).
    void SetMetricDescrip(std::string const& m) { metricDescription_ = m; }
  private:
    Cframes frameToIdx_; ///< Hold indices for all cached frames, -1 for non-cached.
    std::string metricDescription_; ///< Optional description of metric used to gen matrix.
    int sieve_;                     ///< Optional sieve value used to select frames.
};
#endif

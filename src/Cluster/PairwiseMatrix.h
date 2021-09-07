#ifndef INC_CLUSTER_PAIRWISE_MATRIX_H
#define INC_CLUSTER_PAIRWISE_MATRIX_H
//#inc lude "../DataSet_PairwiseCache.h"
//#inc lude "Metric.h"
class DataSet_PairwiseCache;
namespace Cpptraj {
namespace Cluster {
class Cframes;
class Metric;
/// Used to calculate/caching pairwise distances according to a given metric.
class PairwiseMatrix {
  public:
    PairwiseMatrix() : cache_(0), metric_(0) {}

    /// Set up with given metric and optional cache.
    int Setup(Metric*, DataSet_PairwiseCache*);

    // -------------------------------------------
    /// \return distance between given frames.
    double Frame_Distance(int, int) const;
    /// Request that distances for the specified frames be cached.
    int CacheDistances(Cframes const&, int, bool);

    // -------------------------------------------
    //bool HasMetric()           const { return (metric_ != 0); }
    /// \return internal metric, const.
    Metric const& DistMetric() const { return *metric_;       }
    /// \return Pointer to distance metric
    Metric* MetricPtr() const { return metric_; }

    /// \return true if PairwiseMatrix contains a cache.
    bool HasCache()                      const { return (cache_ != 0); }
    /// \return internal cache, const.
    DataSet_PairwiseCache const& Cache() const { return *cache_; }
  private:
    /// Internal routine used to cache pairwise distances.
    int CalcFrameDistances(Cframes const&);

    DataSet_PairwiseCache* cache_;  ///< Hold any cached pairwise distances. 
    Metric* metric_;                ///< The current distance metric.
};

} /* END namespace Cluster */
} /* END namespace Cpptraj */
#endif

#ifndef INC_CLUSTER_PAIRWISE_MATRIX_H
#define INC_CLUSTER_PAIRWISE_MATRIX_H
#include "../DataSet_PairwiseCache.h"
#include "Metric.h"
namespace Cpptraj {
namespace Cluster {

/// Used to calculate/caching pairwise distances according to a given metric.
class PairwiseMatrix {
  public:
    PairwiseMatrix() : cache_(0), metric_(0) {}

    /// Set up with given metric and optional cache.
    int Setup(Metric*, DataSet_PairwiseCache*);

    // -------------------------------------------
    /// \return distance between given frames.TODO const?
    double Frame_Distance(int, int) const;
    /// Request that distances for the specified frames be cached.
    int CacheDistances(Cframes const&);
    /// Print only cached distances. TODO const?
    //virtual void PrintCached() const = 0;

    // -------------------------------------------
    //bool HasMetric()           const { return (metric_ != 0); }
    /// \return internal metric, const.
    Metric const& DistMetric() const { return *metric_;       }
    /// \return internal metric.
//    Metric&       DistMetric()       { return *metric_; }
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

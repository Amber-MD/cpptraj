#ifndef INC_CLUSTER_PAIRWISE_MATRIX_H
#define INC_CLUSTER_PAIRWISE_MATRIX_H
#include "Metric.h"
namespace Cpptraj {
namespace Cluster {

/// Interface for calculating/caching pairwise distances according to a given metric.
class PairwiseMatrix {
  public:
    PairwiseMatrix() : metric_(0) {}
    virtual ~PairwiseMatrix() {}
    // -------------------------------------------
    /// \return distance between given cached frames.
    virtual double GetFdist(int, int) const = 0;
    /// \return distance between given frames.
    virtual double Frame_Distance(int, int) const = 0;
    /// Request that distances for the specified frames be cached.
    virtual int CacheDistances(Cframes const&) = 0;
    // -------------------------------------------
    /// \return internal metric, const.
    Metric const& DistMetric() const { return *metric_; }
    /// \return internal metric.
    Metric&       DistMetric()       { return *metric_; }
  protected:
    /// Used to cache distances; expect internal indices, not absolute cluster frames.
    virtual void SetElement(int, int, double) = 0;
    // -------------------------------------------
    /// Internal routine used to cache pairwise distances.
    int CalcFrameDistances(Cframes const&);
    /// \return Pointer to distance metric
    Metric* distMetric() const { return metric_; }
  private:
    Metric* metric_; ///< The current distance metric.
};

} /* END namespace Cluster */
} /* END namespace Cpptraj */
#endif

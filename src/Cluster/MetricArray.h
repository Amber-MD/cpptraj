#ifndef CPPTRAJ_CLUSTER_METRICARRAY_H
#define CPPTRAJ_CLUSTER_METRICARRAY_H
#include <vector>
namespace Cpptraj {
namespace Cluster {
class Metric;
/// Hold Metrics of various types
class MetricArray {
  public:
    /// CONSTRUCTOR
    MetricArray();
    /// COPY
    MetricArray(MetricArray const&);
    /// ASSIGN
    MetricArray& operator=(MetricArray const&);
    /// Types of distance calculation
    enum DistanceType { MANHATTAN = 0, EUCLID };
  private:
    std::vector<Metric*> metrics_; ///< Hold each Metric
    DistanceType type_;            ///< Type of distance calc to perform
};

}
}
#endif

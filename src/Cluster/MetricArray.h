#ifndef CPPTRAJ_CLUSTER_METRICARRAY_H
#define CPPTRAJ_CLUSTER_METRICARRAY_H
#include <vector>
#include "Metric.h" // Metric::Type
class ArgList;
class DataSet;
class DataSetList;
namespace Cpptraj {
namespace Cluster {

/// Hold Metrics of various types
class MetricArray {
  public:
    /// CONSTRUCTOR
    MetricArray();
    /// DESTRUCTOR
    ~MetricArray();
    /// COPY
    MetricArray(MetricArray const&);
    /// ASSIGN
    MetricArray& operator=(MetricArray const&);
    /// Types of distance calculation
    enum DistanceType { MANHATTAN = 0, EUCLID };
    /// Initialize with data sets and user arguments
    int InitMetricArray(DataSetList const&, ArgList&);
  private:
    /// Clear array
    void Clear();
    /// Allocate metric of given type
    Metric* AllocateMetric(Metric::Type);

    std::vector<Metric*> metrics_; ///< Hold each Metric
    std::vector<DataSet*> sets_;   ///< Sets corresponding to each Metric
    DistanceType type_;            ///< Type of distance calc to perform
};

}
}
#endif

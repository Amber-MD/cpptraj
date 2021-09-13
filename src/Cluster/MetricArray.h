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
    /// Recognized keywords
    static const char* MetricArgs_;
    /// Initialize with data sets and user arguments
    int InitMetricArray(DataSetList const&, ArgList&, int);
    /// Set up all metrics
    int Setup();

    /// \return True if no metrics have been initialized
    bool empty() const { return metrics_.empty(); }
    /// Call Info() routines for all metrics
    void Info() const;
    /// \return Number of points covered by each Metric
    unsigned int Ntotal() const { return ntotal_; }
  private:
    /// Clear array
    void Clear();
    /// Allocate metric of given type
    Metric* AllocateMetric(Metric::Type);

    std::vector<Metric*> metrics_; ///< Hold each Metric
    std::vector<DataSet*> sets_;   ///< Sets corresponding to each Metric
    std::vector<double> weights_;  ///< Weight of each metric
    DistanceType type_;            ///< Type of distance calc to perform
    unsigned int ntotal_;          ///< Total number of points covered by any Metric
};

}
}
#endif

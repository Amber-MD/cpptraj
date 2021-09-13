#ifndef CPPTRAJ_CLUSTER_METRICARRAY_H
#define CPPTRAJ_CLUSTER_METRICARRAY_H
#include <vector>
#include "Metric.h" // Metric::Type
class ArgList;
class DataSet;
class DataSetList;
namespace Cpptraj {
namespace Cluster {
class CentroidArray;
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

    /// Calculate new centroid for each metric
    void NewCentroid(CentroidArray&, Cframes const&);
    /// Calculate centroids for each metric
    void CalculateCentroid(CentroidArray&, Cframes const&);
    /// Update centroids by performing given operation between given frame and centroids.
    void FrameOpCentroid(int, CentroidArray&, double, Metric::CentOpType);
    /// \return distance between given frame and centroids
    double FrameCentroidDist(int, CentroidArray const&);
  private:
    /// Clear array
    void Clear();
    /// Allocate metric of given type
    Metric* AllocateMetric(Metric::Type);
    /// Manhattan distance
    double Dist_Manhattan() const;
    /// Euclidean distance
    double Dist_Euclidean() const;
    /// Distance based on type_
    double DistCalc() const;

    std::vector<Metric*> metrics_; ///< Hold each Metric TODO deal with OpenMP
    std::vector<DataSet*> sets_;   ///< Sets corresponding to each Metric
    std::vector<double> weights_;  ///< Weight of each metric
    std::vector<double> temp_;     ///< For calculations; hold distances from each metric.
    DistanceType type_;            ///< Type of distance calc to perform
    unsigned int ntotal_;          ///< Total number of points covered by any Metric
};

}
}
#endif

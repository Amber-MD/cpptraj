#ifndef CPPTRAJ_CLUSTER_METRICARRAY_H
#define CPPTRAJ_CLUSTER_METRICARRAY_H
#include <vector>
#include "Metric.h" // Metric::Type
class ArgList;
class DataFileList;
class DataSet;
class DataSet_PairwiseCache;
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
    /// Recognized metric keywords
    static const char* MetricArgs_;
    /// Recognized pairwise keywords 1
    static const char* PairwiseArgs1_;
    /// Recognized pairwise keywords 2
    static const char* PairwiseArgs2_;

    /// Initialize with sets to cluster, master data set/file lists (for cache), and user arguments
    int Initialize(DataSetList const&, DataSetList&, DataFileList&, ArgList&, int);
    /// Set up all metrics
    int Setup();
    /// Request that the pairwise cache be filled with distances for specified frames.
    int CacheDistances(Cframes const&, int);

    /// \return True if no metrics have been initialized
    bool empty() const { return metrics_.empty(); }
    /// Call Info() routines for all metrics
    void Info() const;
    /// \return Number of points covered by each Metric
    unsigned int Ntotal() const { return ntotal_; }
    /// \return First Metric (if any) having to do with a COORDS set
    Metric const* CoordsMetric() const;

    // TODO: The CachePtr() and CacheWasAllocated() routines are only needed
    //       because pytraj expects the cluster # vs time set to be
    //       allocated *before* the cache. So it is potentially
    //       allocated via InitMetricArray(), but added to the master
    //       DataSetList later. If/when pytraj is modified to change
    //       this dependency, these two routines can be retired.
    /// \return Pointer to internal pairwise cache.
    DataSet_PairwiseCache* CachePtr() const { return cache_; }
    /// \return true if cache was allocated by MetricArray
    bool CacheWasAllocated()          const { return cacheWasAllocated_; }

    /// Calculate new centroid for each metric
    void NewCentroid(CentroidArray&, Cframes const&);
    /// Calculate centroids for each metric
    void CalculateCentroid(CentroidArray&, Cframes const&);
    /// Update centroids by performing given operation between given frame and centroids.
    void FrameOpCentroid(int, CentroidArray&, double, Metric::CentOpType);
    /// \return distance between given frame and centroids
    double FrameCentroidDist(int, CentroidArray const&);

    /// \return distance between frames (uncached)
    double Uncached_Frame_Distance(int, int);
    /// \return distance between frames (cached or uncached)
    double Frame_Distance(int, int);

  private:
    /// Clear array
    void Clear();
    /// Set up pairwise cache
    int setupPairwiseCache(ArgList&, DataSetList&, DataFileList&);
    /// Allocate metric of given type
    Metric* AllocateMetric(Metric::Type);
    /// Initialize with data sets and user arguments
    int initMetricArray(DataSetList const&, ArgList&);
    /// Manhattan distance
    double Dist_Manhattan(const double*, unsigned int) const;
    /// Euclidean distance
    double Dist_Euclidean(const double*, unsigned int) const;
    /// Distance based on type_
    double DistCalc(std::vector<double> const&) const;
    /// Distance based on type_
    double DistCalc(const double*, unsigned int) const;

    /// \return a string containing description of all metrics
    std::string descriptionOfMetrics() const;
    /// Fill the pairwise cache with distances for specified frames.
    int calcFrameDistances(Cframes const&);

    std::vector<Metric*> metrics_;  ///< Hold each Metric TODO deal with OpenMP
    std::vector<DataSet*> sets_;    ///< Sets corresponding to each Metric
    std::vector<double> weights_;   ///< Weight of each metric
    std::vector<double> temp_;      ///< For calculations; hold distances from each metric.
    int debug_;                     ///< Debug level
    DistanceType type_;             ///< Type of distance calc to perform
    unsigned int ntotal_;           ///< Total number of points covered by any Metric
    DataSet_PairwiseCache* cache_;  ///< Optional cache for frame-frame distances.
    bool cacheWasAllocated_;        ///< True is cache was allocated by InitMetricArray()
    bool pw_mismatch_fatal_;        ///< Controls if PW distances should be recalculated on mismatch
};

}
}
#endif

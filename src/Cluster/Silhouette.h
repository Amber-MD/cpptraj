#ifndef INC_CLUSTER_SILHOUETTE_H
#define INC_CLUSTER_SILHOUETTE_H
#include <vector>
#include <utility> // std::pair
namespace Cpptraj {
namespace Cluster {
class List;
class MetricArray;
/// Cluster silhouette calculation
class Silhouette {
    /// Used to pair frame numbers with silhouette values.
    typedef std::pair<int,double> SilPair;
    /// Used to hold list of frame numbers/silhouette values.
    typedef std::vector<SilPair> SilPairArray;

  public:
    /// CONSTRUCTOR - set debug level
    Silhouette(int);

    /// Used to hold SilPairArray for each cluster
    typedef std::vector<SilPairArray> SilFrameArray;
    /// Generic double array 
    typedef std::vector<double> Darray;

    /// \return Array containing frame silhouette values of each frame for every cluster
    SilFrameArray const& SilhouetteFrameArray() const { return clusterFrameSil_; }
    /// \return Array containing average silhouette values for every cluster
    Darray const& AvgSilhouetteArray() const { return clusterAvgSil_; }

    /// Calculate silhouette values for each cluster
    int CalcSilhouette(List const&, MetricArray&, std::vector<bool> const&, bool);
  private:
    /// For each cluster, hold silhouette values for each frame in the cluster.
    SilFrameArray clusterFrameSil_;
    /// For each cluster, hold avg. silhouette value for cluster (TODO s.d. as well?)
    Darray clusterAvgSil_;
    /// debug level
    int debug_;
};
}
}
#endif

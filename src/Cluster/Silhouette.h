#ifndef INC_CLUSTER_SILHOUETTE_H
#define INC_CLUSTER_SILHOUETTE_H
#include <vector>
#include <utility> // std::pair
#include "../OnlineVarT.h"
class CpptrajFile;
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

    /// Used to determine what indices will be printed out for frame silhouette values
    enum IdxType { IDX_SORTED = 0, IDX_FRAME, IDX_NOT_SPECIFIED };

    /// Used to hold SilPairArray for each cluster
    typedef std::vector<SilPairArray> SilFrameArray;
    /// Used to hold average silhouette value for each cluster 
    typedef std::vector<Stats<double> > Darray;

    /// Initialize with frame silhouette index type
    int Init(IdxType);

    /// \return Array containing frame silhouette values of each frame for every cluster
    SilFrameArray const& SilhouetteFrameArray() const { return clusterFrameSil_; }
    /// \return Array containing average silhouette values for every cluster
    Darray const& AvgSilhouetteArray() const { return clusterAvgSil_; }

    /// Print Silhouette frame values to file.
    int PrintSilhouetteFrames(CpptrajFile&, List const&) const;
    /// Print average Silhouette values to file.
    int PrintAvgSilhouettes(CpptrajFile&, List const&) const;

    /// Calculate silhouette values for each cluster
    int CalcSilhouette(List const&, MetricArray&, std::vector<bool> const&, bool);
  private:
    /// \return 1 if given # clusters does not match what is in Silhouette
    int numMismatchErr(const char*, unsigned int) const;

    /// For sorting cluster frame silhouettes by silhouette value.
    struct sort_by_sil_val {
      inline bool operator()(SilPair const& p0, SilPair const& p1)
      {
        if (p0.second == p1.second)
          return (p0.first < p1.first);
        else
          return (p0.second < p1.second);
      }
    };

    /// For each cluster, hold silhouette values for each frame in the cluster.
    SilFrameArray clusterFrameSil_;
    /// For each cluster, hold avg. silhouette value for cluster (TODO s.d. as well?)
    Darray clusterAvgSil_;
    /// debug level
    int debug_;
    /// Control what indices will be used when printing frame silhouette values
    IdxType silIdxType_;
};
}
}
#endif

#ifndef INC_CLUSTER_ALGORITHM_KMEANS_H
#define INC_CLUSTER_ALGORITHM_KMEANS_H
#include "Algorithm.h"
#include "../Random.h"
namespace Cpptraj {
namespace Cluster {

/// Implement K-means clustering
class Algorithm_Kmeans : public Algorithm {
  public:
    Algorithm_Kmeans();
    static void Help();
    int Setup(ArgList&);
    void Info() const;
    void Results(CpptrajFile&) const;
    int DoClustering(List&, Cframes const&, MetricArray&);
    void Timing(double) const {}
  private:
    typedef std::vector<int> Iarray;
    enum KmeansModeType { SEQUENTIAL, RANDOM };
    enum KmeansInitType { INIT_SPREAD=0, INIT_KMPP };

    Iarray FindSeedsFromClusters(List&, MetricArray&) const;
    int FindKmeansSeeds(Cframes const&, MetricArray&);
    /// Kmeans++ initialization
    int init_kmeanspp(List&, Cframes const&, MetricArray&);

    Random_Number RN_;
    int nclusters_; ///< Target number of clusters.
    int kseed_;
    int maxIt_;
    Iarray SeedIndices_; ///< Hold indices into framesToCluster of seed frames
    KmeansModeType mode_;
    bool clusterToClusterCentroid_;
    KmeansInitType init_; ///< Initialization strategy
};

}
}
#endif

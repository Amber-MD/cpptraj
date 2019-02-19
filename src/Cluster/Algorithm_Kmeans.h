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
    int DoClustering(List&, Cframes const&, PairwiseMatrix const&);
    void Timing(double) const {}
  private:
    typedef std::vector<int> Iarray;
    enum KmeansModeType { SEQUENTIAL, RANDOM };

    int FindKmeansSeeds(Cframes const&, PairwiseMatrix const&);
    void ShufflePoints(Iarray&);

    Random_Number RN_;
    int nclusters_; ///< Target number of clusters.
    int kseed_;
    int maxIt_;
    Iarray SeedIndices_;
    KmeansModeType mode_;
    bool clusterToClusterCentroid_;
};

}
}
#endif

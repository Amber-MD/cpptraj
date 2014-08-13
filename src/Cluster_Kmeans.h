#ifndef INC_CLUSTER_KMEANS_H
#define INC_CLUSTER_KMEANS_H
#include "ClusterList.h"
class Cluster_Kmeans : public ClusterList {
  public:
    Cluster_Kmeans();
    static void Help();
    int SetupCluster(ArgList&);
    void ClusteringInfo();
    int Cluster();
    void AddSievedFrames();
    void ClusterResults(CpptrajFile&) const;
  private:
    enum KmeansModeType { SEQUENTIAL, RANDOM };

    int FindKmeansSeeds();

    int nclusters_; ///< Target number of clusters.
    int kseed_;
    int maxIt_;
    typedef std::vector<int> Iarray;
    Iarray SeedIndices_;
    Iarray FramesToCluster_;
    KmeansModeType mode_;
};
#endif

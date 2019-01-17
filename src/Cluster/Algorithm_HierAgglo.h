#ifndef INC_CLUSTER_ALGORITHM_HIERAGGLO_H
#define INC_CLUSTER_ALGORITHM_HIERAGGLO_H
#include "Algorithm.h"
#include "../Timer.h"
#include "../ClusterMatrix.h"
namespace Cpptraj {
namespace Cluster {

/// Implement hierarchical agglomerative clustering
class Algorithm_HierAgglo : public Algorithm {
  public:
    Algorithm_HierAgglo();
    static void Help();
    int Setup(ArgList&);
    void Info() const;
    void Results(CpptrajFile&) const;
    int DoClustering(List&, Cframes const&, PairwiseMatrix const&);
    void Timing(double) const;
  private:
    void buildInitialClusters(List&, Cframes const&, Metric&);
    void InitializeClusterDistances();
    int MergeClosest(List&);
    void calcMinDist(List::cluster_it&);
    void calcMaxDist(List::cluster_it&);
    void calcAvgDist(List::cluster_it&);

    /// Type of distance calculation between clusters.
    enum LINKAGETYPE  { SINGLELINK = 0, AVERAGELINK, COMPLETELINK };
    int nclusters_;       ///< Target # of clusters.
    double epsilon_;      ///< Once the min distance between clusters is > epsilon, stop.
    LINKAGETYPE linkage_; ///< Cluster Linkage type.
    CpptrajFile eps_v_n_; ///< Write epsilon vs # clusters.
    ClusterMatrix ClusterDistances_;
#   ifdef TIMER
    Timer time_findMin_;
    Timer time_mergeFrames_;
    Timer time_calcLinkage_;
#   endif
};

} /* END namespace Cluster */
} /* END namespace Cpptraj */
#endif

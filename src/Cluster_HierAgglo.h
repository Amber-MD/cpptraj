#ifndef INC_CLUSTER_HIERAGGLO_H
#define INC_CLUSTER_HIERAGGLO_H
#include "ClusterList.h"
#include "ClusterMatrix.h"
#ifdef TIMER
# include "Timer.h"
#endif
class Cluster_HierAgglo : public ClusterList {
  public:
    Cluster_HierAgglo();
    static void Help();
    int SetupCluster(ArgList&);
    void ClusteringInfo() const;
    int Cluster();
    /// \return Distance between given clusters based on current linkage
    double ClusterDistance(ClusterNode const&, ClusterNode const&) const;
#   ifdef TIMER
    void Timing(double) const;
#   endif
    void AddSievedFrames() { AddSievedFramesByCentroid(); }
    void ClusterResults(CpptrajFile&) const;
  private:
    void InitializeClusterDistances();
    int MergeClosest();
    void calcMinDist(cluster_it&);
    void calcMaxDist(cluster_it&);
    void calcAvgDist(cluster_it&);

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
#endif

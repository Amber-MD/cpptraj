#ifndef INC_CLUSTER_HIERAGGLO_H
#define INC_CLUSTER_HIERAGGLO_H
#include "ClusterList.h"
#ifdef TIMER
# include "Timer.h"
#endif
class Cluster_HierAgglo : public ClusterList {
  public:
    /// Type of distance calculation between clusters.
    enum LINKAGETYPE  { SINGLELINK = 0, AVERAGELINK, COMPLETELINK };
    Cluster_HierAgglo();
    static void Help();
    int SetupCluster(ArgList&);
    void ClusteringInfo() const;
    int Cluster();
#   ifdef TIMER
    void Timing(double) const;
#   endif
    void AddSievedFrames() { AddSievedFramesByCentroid(); }
    void ClusterResults(CpptrajFile&) const;
  private:
    int nclusters_;       ///< Target # of clusters.
    double epsilon_;      ///< Once the min distance between clusters is > epsilon, stop.
    LINKAGETYPE linkage_; ///< Cluster Linkage type.
    CpptrajFile eps_v_n_; ///< Write epsilon vs # clusters.
#   ifdef NEWCODE
    /** Upper-triangle matrix containing sum of distances between frames in
      * cluster i to frames in cluster j, for use with average linkage. */
    Matrix<double> SumDistToCluster_;
#   endif
#   ifdef TIMER
    Timer time_findMin_;
    Timer time_mergeFrames_;
    Timer time_calcLinkage_;
#   endif
    void InitializeClusterDistances();
    int MergeClosest();
    void calcMinDist(cluster_it&);
    void calcMaxDist(cluster_it&);
    void calcAvgDist(cluster_it&);
};
#endif

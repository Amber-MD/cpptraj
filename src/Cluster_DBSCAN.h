#ifndef INC_CLUSTER_DBSCAN_H
#define INC_CLUSTER_DBSCAN_H
#include "ClusterList.h"
class Cluster_DBSCAN : public ClusterList {
  public:
    Cluster_DBSCAN();
    static void Help();
    int SetupCluster(ArgList&);
    void ClusteringInfo();
    int Cluster();
  private:
    int minPoints_;  ///< Min # of points needed to make a cluster.
    double epsilon_; ///< Distance criterion for cluster formation.

    void RegionQuery(std::vector<int>&, std::vector<int> const&, int);
};
#endif

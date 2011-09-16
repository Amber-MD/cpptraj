#ifndef INC_CLUSTERLIST_H
#define INC CLUSTERLIST_H
/// Class: ClusterList
/// Store information on clusters
#include <list>
#include "TriangleMatrix.h"
#include "DataSet.h"
class ClusterList {
  public :
    // Store individual cluster info; frame numbers, centroid, etc.
    struct clusterNode {
      std::list<int> frameList; // List of frames in the cluster
      int num;                  // Cluster number (index in Distances)
    };
    std::list<clusterNode> clusters;
    int maxframes; // Total number of frames being clustered

    struct cluster_cmp {
      bool operator()(clusterNode first, clusterNode second) const {
        if (first.frameList.size() > second.frameList.size())
          return true;
        else
          return false;
      }
    };

    ClusterList();
    ~ClusterList();
   
    int AddCluster(std::list<int> *, int);
    void PrintClusters();
    std::list<clusterNode>::iterator GetCluster(int);
    int Merge(std::list<clusterNode>::iterator,std::list<clusterNode>::iterator);
    void calcMinDist(std::list<ClusterList::clusterNode>::iterator,
                     TriangleMatrix *, TriangleMatrix *);
    void calcAvgDist(std::list<ClusterList::clusterNode>::iterator,
                     TriangleMatrix *, TriangleMatrix *);
    void Cnumvtime(DataSet *);
    void Renumber();
    void Summary(char *);
};
#endif

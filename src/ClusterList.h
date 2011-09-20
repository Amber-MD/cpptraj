#ifndef INC_CLUSTERLIST_H
#define INC CLUSTERLIST_H
/// Class: ClusterList
/// This class holds all the individual clusters, as well as routines that
/// can be used to perform clustering (metric recalculation, cluster merging,
/// and so on). The distance calculation routines require that a
/// triangle matrix with distances between all frames be previously
/// calculated.
#include <list>
#include "TriangleMatrix.h"
class ClusterList {
  public:
    enum LINKAGETYPE {SINGLELINK, AVERAGELINK, COMPLETELINK};

  private:
    int debug;
    // clusterNode: Store individual cluster info; frame numbers, centroid, etc.
    struct clusterNode {
      std::list<int> frameList; // List of frames in the cluster
      int num;                  // Cluster number (index in Distances)
      int centroid;             // Frame number of the cluster centroid
    };
    std::list<clusterNode> clusters;

    int maxframes;                   // Total number of frames being clustered
    TriangleMatrix *FrameDistances;  // Distances between each frame
    TriangleMatrix ClusterDistances; // Distances between each cluster
    LINKAGETYPE Linkage;             // Type of distance calculation for clusters 

    // cluster_cmp: Used to sort the cluster list
    struct cluster_cmp {
      bool operator()(clusterNode first, clusterNode second) const {
        if (first.frameList.size() > second.frameList.size())
          return true;
        else
          return false;
      }
    };

    int Merge(std::list<clusterNode>::iterator,std::list<clusterNode>::iterator);
    void FindCentroid(std::list<clusterNode>::iterator);
    std::list<clusterNode>::iterator GetClusterIt(int);

    // Distance calculation routines
    void calcMinDist(std::list<clusterNode>::iterator);
    void calcMaxDist(std::list<clusterNode>::iterator);
    void calcAvgDist(std::list<clusterNode>::iterator);

    std::list<ClusterList::clusterNode>::iterator currentCluster;

  public :
    ClusterList();
    ~ClusterList();
  
    void SetLinkage(LINKAGETYPE Lin) { Linkage = Lin; }
    int Nclusters() { return (int) clusters.size(); }

    void Initialize(TriangleMatrix *);
    int AddCluster(std::list<int> *, int);
    void PrintClusters();
    int MergeClosest(double);
    void Renumber();
    void Summary(char *);

    void Begin();
    bool End();
    void NextCluster();
    int CurrentNum();
    int CurrentCentroid();
    std::list<int>::iterator CurrentFrameBegin();
    std::list<int>::iterator CurrentFrameEnd();
};
#endif

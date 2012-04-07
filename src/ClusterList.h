#ifndef INC_CLUSTERLIST_H
#define INC CLUSTERLIST_H
#include "ClusterNode.h" 
#include "TriangleMatrix.h"
// Class: ClusterList
/** This class holds all the individual clusters, as well as routines that
  * can be used to perform clustering (metric recalculation, cluster merging,
  * and so on). The distance calculation routines require that a
  * triangle matrix with distances between all frames be previously
  * calculated.
  */
class ClusterList {
  public:
    enum LINKAGETYPE {SINGLELINK, AVERAGELINK, COMPLETELINK};

    ClusterList();
 
    void SetLinkage(LINKAGETYPE Lin) { Linkage_ = Lin; }
    int Nclusters() { return (int)clusters_.size(); }

    void SetDebug(int);
    void Initialize(TriangleMatrix *);
    int AddCluster(std::list<int>&, int);
    void PrintClusters();
    void PrintClustersToFile(char *);
    void PrintRepFrames();
    int MergeClosest(double);
    void Renumber();
    void Summary(char *);
    void Summary_Half(char *);
    bool CheckEpsilon(double);

/*    void Begin();
    bool End();
    void NextCluster();
    int CurrentNum();
    int CurrentCentroid();
    std::list<int>::iterator CurrentFrameBegin();
    std::list<int>::iterator CurrentFrameEnd();*/
  private:
    static const char XMGRACE_COLOR[][12];

    int debug_;
    // clusterNode: Store individual cluster info; frame numbers, centroid, etc.
    std::list<ClusterNode> clusters_;
    /// Total number of frames being clustered
    int maxframes_;
    /// Distances between each frame
    TriangleMatrix *FrameDistances_;
    /// Distances between each cluster 
    TriangleMatrix ClusterDistances_;
    /// Type of distance calculation for clusters 
    LINKAGETYPE Linkage_;

    typedef std::list<ClusterNode>::iterator cluster_it;

    int Merge(cluster_it&, cluster_it&);
    void FindCentroid(cluster_it&);
    //std::list<clusterNode>::iterator GetClusterIt(int);

    // Distance calculation routines
    void calcMinDist(cluster_it&);
    void calcMaxDist(cluster_it&);
    void calcAvgDist(cluster_it&);

    void CalcEccentricity(cluster_it&);
};
#endif

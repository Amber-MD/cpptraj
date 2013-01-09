#ifndef INC_CLUSTERLIST_H
#define INC CLUSTERLIST_H
#include "ClusterNode.h" 
// Class: ClusterList
/** This class holds all the individual clusters, as well as routines that
  * can be used to perform clustering (metric recalculation, cluster merging,
  * and so on).
  */
class ClusterList {
  public:
    /// Type of distance calculation between clusters.
    enum LINKAGETYPE  { SINGLELINK = 0, AVERAGELINK, COMPLETELINK };
    enum DistModeType { USE_FRAMES = 0, USE_FILE };
    ClusterList();
    ~ClusterList();
    int Nclusters()                  const { return (int)clusters_.size(); }

    void SetDebug(int);
    void Renumber();
    void Summary(std::string const&,int);
    void Summary_Half(std::string const&,int);
    void PrintClustersToFile(std::string const&,int);
    void PrintClusters();

    int CalcFrameDistances(std::string const&, DataSet*, DistModeType, 
                           bool, bool, bool, std::string const&, int);
    void AddSievedFrames();
    // Clustering Methods
    int ClusterHierAgglo(double, int, LINKAGETYPE);

    // Const Iterator over clusters
    typedef std::list<ClusterNode>::const_iterator cluster_iterator;
    cluster_iterator begincluster() { return clusters_.begin(); }
    cluster_iterator endcluster()   { return clusters_.end();   }
  private:
    /// Iterator over clusters
    typedef std::list<ClusterNode>::iterator cluster_it;
    static const char* XMGRACE_COLOR[];
    int debug_;
    /// Store individual cluster info; frame numbers, centroid, etc.
    std::list<ClusterNode> clusters_;
    /// Distances between each frame.
    ClusterMatrix FrameDistances_;
    /// Distances between each cluster.
    ClusterMatrix ClusterDistances_;
    /// Used to calculate distances between frames and/or centroids.
    ClusterDist* Cdist_;

    int AddCluster(std::list<int> const&);

    // Hierarchical Agglomerative calculation routines
    void InitializeClusterDistances(LINKAGETYPE);
    int MergeClosest(double, LINKAGETYPE);
    void calcMinDist(cluster_it&);
    void calcMaxDist(cluster_it&);
    void calcAvgDist(cluster_it&);

    double ComputeDBI(CpptrajFile&);
};
#endif

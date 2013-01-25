#ifndef INC_CLUSTERLIST_H
#define INC_CLUSTERLIST_H
#include "ArgList.h"
#include "ClusterNode.h" 
// Class: ClusterList
/** This base class holds all the individual clusters, as well as routines 
  * that can be used to obtain information on clusters after clustering.
  */
class ClusterList {
  public:
    enum DistModeType { USE_FRAMES = 0, USE_FILE };
    ClusterList();
    virtual ~ClusterList();
    int Nclusters()                  const { return (int)clusters_.size(); }

    void SetDebug(int);
    void Renumber();
    void Summary(std::string const&,int);
    void Summary_Half(std::string const&,int,int);
    void PrintClustersToFile(std::string const&,int);
    void PrintClusters();

    int CalcFrameDistances(std::string const&, ClusterDist::DsArray const&, DistModeType, 
                           bool, bool, bool, std::string const&, int);
    void AddSievedFrames();
    // Inherited by individual clustering methods
    virtual int SetupCluster(ArgList&) = 0;
    virtual void ClusteringInfo() = 0;
    virtual int Cluster() = 0;

    // Const Iterator over clusters
    typedef std::list<ClusterNode>::const_iterator cluster_iterator;
    const cluster_iterator begincluster() const { return clusters_.begin(); }
    const cluster_iterator endcluster()   const { return clusters_.end();   }
  protected:
    /// Iterator over clusters
    typedef std::list<ClusterNode>::iterator cluster_it;
    int debug_;
    /// Store individual cluster info; frame numbers, centroid, etc.
    std::list<ClusterNode> clusters_;
    /// Distances between each frame.
    ClusterMatrix FrameDistances_;
    /// Distances between each cluster.
    ClusterMatrix ClusterDistances_;
    /// Used to calculate distances between frames and/or centroids.
    ClusterDist* Cdist_;
    /// Add specified frames to a new cluster.
    int AddCluster(ClusterDist::Cframes const&);
    /// Calculate the Davies-Bouldin index of clusters.
    double ComputeDBI(CpptrajFile&);
  private:
    static const char* XMGRACE_COLOR[];
};
#endif

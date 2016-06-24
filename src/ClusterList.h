#ifndef INC_CLUSTERLIST_H
#define INC_CLUSTERLIST_H
#include <list>
#include "ArgList.h"
#include "ClusterNode.h"
// Class: ClusterList
/** This base class holds all the individual clusters, as well as routines 
  * that can be used to obtain information on clusters after clustering.
  */
class ClusterList {
  public:
    enum DistMetricType { RMS = 0, DME, SRMSD, DATA };
    static const char* MetricString( DistMetricType );
    ClusterList();
    virtual ~ClusterList();
    int Nclusters()                  const { return (int)clusters_.size(); }

    void SetDebug(int);
    /// Add back sieved frames, update centroids, sort by cluster population.
    void Renumber(bool);
    /// Determine which frames in each cluster are best representative using cumulative distance.
    int FindBestRepFrames_CumulativeDist();
    /// Determine which frames in each cluster are best representative by distance to centroid.
    int FindBestRepFrames_Centroid();
    /// Print overall summary of clusters.
    void Summary(std::string const&,int) const;
    /// Print summary of clusters separated by parts.
    void Summary_Part(std::string const&,int,std::vector<int> const&) const;
    /// Print cluster info.
    void PrintClustersToFile(std::string const&,int) const;
    /// DEBUG: Print clusters to STDOUT
    void PrintClusters() const;
    /// Set up appropriate cluster distance calculation
    int SetupCdist( ClusterDist::DsArray const&, DistMetricType, bool, bool, std::string const&);
    /// Calculate distances between frames if necessary.
    int CalcFrameDistances(DataSet*, ClusterDist::DsArray const&, int, int);
    // ----- Inherited by individual clustering methods ----
    /// Process algorithm-specific keywords
    virtual int SetupCluster(ArgList&) = 0;
    /// Print summary of clustering to be performed
    virtual void ClusteringInfo() const = 0;
    /// Perform clustering.
    virtual int Cluster() = 0;
    /// \return distance between given clusters. Default is distance between centroids.
    virtual double ClusterDistance(ClusterNode const&, ClusterNode const&) const;
#   ifdef TIMER
    virtual void Timing(double) const = 0;
#   endif
    /// Iterator over clusters
    typedef std::list<ClusterNode>::iterator cluster_it;
    cluster_it begin() { return clusters_.begin(); }
    cluster_it end()   { return clusters_.end();   }
    /// Const Iterator over clusters
    typedef std::list<ClusterNode>::const_iterator cluster_iterator;
    const cluster_iterator begincluster() const { return clusters_.begin(); }
    const cluster_iterator endcluster()   const { return clusters_.end();   }
    /// Remove clusters with no members.
    void RemoveEmptyClusters();
    /// Calculate cluster silhouettes
    void CalcSilhouette(std::string const&) const;

    void DrawGraph(bool,DataSet*,double,int) const;
  protected:
    /// Restore sieved frames to clusters
    virtual void AddSievedFrames() = 0;
    virtual void ClusterResults(CpptrajFile&) const = 0;

    /// \return Distance between specified frames. Use FrameDistances if frames were not sieved.
    inline double Frame_Distance(int,int) const;
    /// Add each sieved frame to the nearest cluster based on frame to centroid distance.
    void AddSievedFramesByCentroid();
    DataSet_Cmatrix const& FrameDistances() const { return *frameDistances_; }
    int debug_;
    /// Store individual cluster info; frame numbers, centroid, etc.
    std::list<ClusterNode> clusters_;
    /// Used to calculate distances between frames and/or centroids.
    ClusterDist* Cdist_;
    /// Add specified frames to a new cluster.
    int AddCluster(ClusterDist::Cframes const&);
  private:
    static const char* XMGRACE_COLOR[];
    /// Determine max name width
    unsigned int DetermineNameWidth() const;
    /// Calculate the Davies-Bouldin index of clusters. Centroids must be up-to-date.
    double ComputeDBI(CpptrajFile&) const;
    /// Calculate pseudo-F statistic.
    double ComputePseudoF(CpptrajFile&) const;

    /// Hold pointer to matrix containing distances between each frame.
    DataSet_Cmatrix* frameDistances_;
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
double ClusterList::Frame_Distance(int f1, int f2) const {
  if (FrameDistances().FrameWasSieved(f1) ||
      FrameDistances().FrameWasSieved(f2))
    return Cdist_->FrameDist(f1, f2);
  else
    return FrameDistances().GetFdist(f1, f2);
}
#endif

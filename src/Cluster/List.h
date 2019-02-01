#ifndef INC_CLUSTER_LIST_H
#define INC_CLUSTER_LIST_H
#include <list>
#include "Node.h"
#include "../DataSet_integer.h"
namespace Cpptraj {
namespace Cluster {

/// Hold all individual clusters.
/** Currently implemented as an STL list since sorting and erasing are more
  * efficient.
  */
class List {
  public:
    List() {}
    /// Set debug level
    void SetDebug(int d) { debug_ = d; }
    /// Iterator over clusters
    typedef std::list<Node>::iterator cluster_it;
    /// Iterator to beginning
    cluster_it begin() { return clusters_.begin(); }
    /// Iterator to end
    cluster_it end()   { return clusters_.end();   }
    /// Const iterator over clusters
    typedef std::list<Node>::const_iterator cluster_iterator;
    /// Const iterator to beginning
    const cluster_iterator begincluster() const { return clusters_.begin(); }
    /// Const iterator to end
    const cluster_iterator endcluster()   const { return clusters_.end();   }
    /// \return first cluster
    Node const& front()                   const { return clusters_.front(); }
    /// \return current number of clusters.
    int Nclusters()        const { return (int)clusters_.size(); }
    /// \return true if no clusters
    bool empty()           const { return clusters_.empty(); }
    /// \return Array containing noise points
    Cframes const& Noise() const { return noise_; }
    /// Print clusters to stdout
    void PrintClusters() const;
    /// Add new cluster
    void AddCluster( Node const& n )     { clusters_.push_back( n ); }
    /// Remove existing cluster via iterator
    void RemoveCluster( cluster_it& it ) { clusters_.erase( it ); }
    /// Generate cluster number vs time data set
    int CreateCnumVsTime(DataSet_integer*, unsigned int) const;
    /// Sort clusters by population and renumber.
    int Sort();
    /// Add given frame as noise.
    void AddNoise(int f) { noise_.push_back( f ); }
    /// Update centroids TODO check if they need updating
    void UpdateCentroids(Metric*);
    /// Add given frames to clusters based on distance to centroid. TODO save original frames
    void AddFramesByCentroid(Cframes const&, Metric*);
    /// Add given frames to clusters based on distance to centroid and cutoff.
    void AddFramesByCentroid(Cframes const&, Metric*, bool, double);
    /// Calculate the Davies-Bouldin index.
    double ComputeDBI(std::vector<double>&, Metric*) const;
    /// Calculate pseudo-F
    double ComputePseudoF(double&, Metric*) const;
    /// Calculate cluster and cluster frame silhouettes
    int CalcSilhouette(PairwiseMatrix const&, Cframes const&, bool);
  private:
    typedef std::list<Node> Narray;
    Narray clusters_; ///< Hold all clusters.
    Cframes noise_;   ///< Hold any frames classified as noise.
    int debug_;
};

} /* END namespace Cluster */
} /* END namespace Cpptraj */

#endif

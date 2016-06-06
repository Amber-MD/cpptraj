#ifndef INC_CLUSTER_DBSCAN_H
#define INC_CLUSTER_DBSCAN_H
#include "ClusterList.h"
/** Ester, Kriegel, Sander, Xu; Proceedings of 2nd International Conference
  * on Knowledge Discovery and Data Mining (KDD-96); pp 226-231.
  */
class Cluster_DBSCAN : public ClusterList {
  public:
    Cluster_DBSCAN();
    static void Help();
    int SetupCluster(ArgList&);
    void ClusteringInfo();
    int Cluster();
    void AddSievedFrames();
    void ClusterResults(CpptrajFile&) const;
  private:
    typedef std::vector<int> Iarray;

    bool ExpandCluster(unsigned int, int);
    void RegionQuery(Iarray&, int) const;
    void ComputeKdist( int, std::vector<int> const& ) const ;
    void ComputeKdistMap( Range const&, std::vector<int> const& ) const ;

    Iarray Status_;        ///< Status of each point: unclassified, noise, or in cluster
    Iarray seeds_;         ///< Results from first RegionQuery
    Iarray result_;        ///< Results from seed RegionQueries
    int minPoints_;        ///< Min # of points needed to make a cluster.
    double epsilon_;       ///< Distance criterion for cluster formation.
    Range kdist_;
    std::string k_prefix_; ///< Kdist output file prefix.
    bool sieveToCentroid_; ///< If true sieve only based on closeness to centroid.
};
#endif

#ifndef INC_CLUSTER_ALGORITHM_DBSCAN_H
#define INC_CLUSTER_ALGORITHM_DBSCAN_H
#include "Algorithm.h"
namespace Cpptraj {
namespace Cluster {

/** Ester, Kriegel, Sander, Xu; Proceedings of 2nd International Conference
  * on Knowledge Discovery and Data Mining (KDD-96); pp 226-231.
  */
class Algorithm_DBscan : public Algorithm {
  public:
    Algorithm_DBscan();
    static void Help();
    int Setup(ArgList&);
    void Info() const;
    void Results(CpptrajFile&) const;
    int DoClustering(List&, Cframes const&, PairwiseMatrix const&);
    void Timing(double) const {}

    //void AddSievedFrames(); // TODO fix this
    double Epsilon() const { return epsilon_; }
  private:
    typedef std::vector<int> Iarray;

    bool ExpandCluster(Cframes const&, PairwiseMatrix const&, unsigned int, int);
    void RegionQuery(Iarray&, Cframes const&, PairwiseMatrix const&, int) const;
    void ComputeKdist( int, Cframes const&, PairwiseMatrix const&) const ;
    void ComputeKdistMap( Range const&, Cframes const&, PairwiseMatrix const& ) const ;

    Iarray Status_;        ///< Status of each point: unclassified, noise, or in cluster
    Iarray seeds_;         ///< Results from first RegionQuery
    Iarray result_;        ///< Results from seed RegionQueries
    int minPoints_;        ///< Min # of points needed to make a cluster.
    double epsilon_;       ///< Distance criterion for cluster formation.
    Range kdist_;
    std::string k_prefix_; ///< Kdist output file prefix.
    bool sieveToCentroid_; ///< If true sieve only based on closeness to centroid.
};

} /* END namespace Cluster */
} /* END namespace Cpptraj */
#endif

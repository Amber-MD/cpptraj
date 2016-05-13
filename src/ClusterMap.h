#ifndef INC_CLUSTERMAP_H
#include "DataSet_2D.h"
#ifdef TIMER
#include "Timer.h"
#endif
class ClusterMap {
  public:
    ClusterMap();
    /// Initialize with given epsilon, min points, and cmap_square
    int Init(double,int,bool);
    /// Perform clustering on given 2D data set.
    int DoCluster(DataSet_2D const&);

    class Cluster;
    typedef std::vector<Cluster> Carray;
    typedef std::vector<int> Iarray;
  private:

    int DoDBSCAN(DataSet_2D const&);
    /// \return true if cluster could be expanded, false if point is noise.
    bool ExpandCluster(unsigned int, int, DataSet_2D const&);
    /// Set given Iarray with neighbors of given point.
    void RegionQuery(Iarray&, int, DataSet_2D const&);
    /// Add given points to a new cluster.
    void AddCluster(Iarray const&, DataSet_2D const&);

    Carray clusters_;  ///< Hold all clusters.
    Iarray Status_;    ///< Status of each point: unclassified, noise, or in cluster
    Iarray seeds_;     ///< Results from first RegionQuery
    Iarray result_;    ///< Results from seed RegionQueries
    double epsilon_;   ///< Distance to search for neighboring points within.
    double epsilon2_;  ///< Epsilon squared.
    double Avg_;       ///< Average value of map, used as cutoff (points below are noise).
    int minPoints_;    ///< Minimum number of points within epsilon to qualify as cluster.
    int idx_offset_;   ///< Max # of rows/cols to offset in RegionQuery() based on epsilon.
    bool cmap_square_; ///< If true write cluster map by min/max rows/cols.
#   ifdef TIMER
    Timer t_overall_;
    Timer t_query1_;
    Timer t_query2_;
#   endif
};
// -----------------------------------------------
class ClusterMap::Cluster {
  public:
    Cluster() : avg_(0.0), cnum_(-1), min_col_(-1), max_col_(-1), min_row_(-1), max_row_(-1) {}
    Cluster(Iarray const& p, double A, int C, int minc, int maxc, int minr, int maxr) :
      points_(p), avg_(A), cnum_(C), min_col_(minc), max_col_(maxc), min_row_(minr), max_row_(maxr)
      {}
    void SetCnum(int c) { cnum_ = c; }
    Iarray const& Points() const { return points_;  }
    double Avg()           const { return avg_;     }
    int Cnum()             const { return cnum_;    }
    int MinCol()           const { return min_col_; }
    int MaxCol()           const { return max_col_; }
    int MinRow()           const { return min_row_; }
    int MaxRow()           const { return max_row_; }
    // Use > since we give higher priority to larger clusters
    bool operator <(Cluster const& rhs) const {
      return ( points_.size() > rhs.points_.size() );
    }
  private:
    Iarray points_;
    double avg_;
    int cnum_;
    int min_col_;
    int max_col_;
    int min_row_;
    int max_row_;
};
#endif

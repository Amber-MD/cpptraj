#ifndef INC_CLUSTER_ALGORITHM_DPEAKS_H
#define INC_CLUSTER_ALGORITHM_DPEAKS_H
#include "Algorithm.h"
#include <string>
#include <vector>
namespace Cpptraj {
namespace Cluster {

/// Implement density peaks clustering algorithm
class Algorithm_DPeaks : public Algorithm {
  public:
    Algorithm_DPeaks();
    static void Help();
    int Setup(ArgList&);
    void Info() const;
    void Results(CpptrajFile&) const;
    int DoClustering(List&, Cframes const&, MetricArray&);
    void Timing(double) const {}
    double Epsilon() const { return epsilon_; }
  private:
    void AssignClusterNum(int, int&);
    int Cluster_GaussianKernel(Cframes const&, MetricArray&);
    int Cluster_DiscreteDensity(Cframes const&, MetricArray&);
    int ChoosePointsAutomatically();
    int ChoosePointsManually();
    int ChoosePointsFromClusters(List const&, Cframes const&);

    enum ChooseType {PLOT_ONLY = 0, MANUAL, AUTOMATIC};
    std::string dvdfile_;
    std::string rafile_;
    std::string radelta_;
    double densityCut_;
    double distanceCut_;
    double epsilon_;
    ChooseType choosePoints_;
    int avg_factor_;
    bool calc_noise_;
    bool useGaussianKernel_;
    // -------------------------------------------
    class Cpoint {
      public:
        Cpoint() :
          dist_(-1.0), density_(0.0), pointsWithinEps_(0), fnum_(-1), 
          nidx_(-1), oidx_(-1), cnum_(-1) {}
        Cpoint(int f, int o) :
          dist_(-1.0), density_(0.0), pointsWithinEps_(0), fnum_(f),
          nidx_(-1), oidx_(o), cnum_(-1) {}
        Cpoint(int f) :
          dist_(-1.0), density_(0.0), pointsWithinEps_(0), fnum_(f),
          nidx_(-1), oidx_(-1), cnum_(-1) {}
        Cpoint(Cpoint const& rhs) :
          dist_(rhs.dist_), density_(rhs.density_),
          pointsWithinEps_(rhs.pointsWithinEps_), fnum_(rhs.fnum_),
          nidx_(rhs.nidx_), oidx_(rhs.oidx_), cnum_(rhs.cnum_) {}
        Cpoint& operator=(Cpoint const& rhs) {
          if (&rhs != this) {
            dist_ = rhs.dist_; density_ = rhs.density_;
            pointsWithinEps_ = rhs.pointsWithinEps_; fnum_ = rhs.fnum_;
            nidx_ = rhs.nidx_; oidx_ = rhs.oidx_; cnum_ = rhs.cnum_;
          }
          return *this;
        }
        /// Used to sort Carray by pointsWithinEpsilon, ascending 
        //bool operator<(Cpoint const& second) const
        struct pointsWithinEps_sort {
          inline bool operator()(Cpoint const& first, Cpoint const& second) const {
            if (first.pointsWithinEps_ == second.pointsWithinEps_)
              return (first.fnum_ < second.fnum_);
            else
              return (first.pointsWithinEps_ < second.pointsWithinEps_);
          }
        };
        /// Sort by density, descending
        struct density_sort_descend {
          inline bool operator()(Cpoint const& first, Cpoint const& second) const {
              return (first.density_ > second.density_);
          }
        };
        /// Used to sort Carray by cluster number
        struct cnum_sort {
          inline bool operator()(Cpoint const& first, Cpoint const& second) const {
            if (first.cnum_ == second.cnum_)
              return (first.fnum_ < second.fnum_);
            else
              return (first.cnum_ < second.cnum_);
          }
        };
        /// Used to sort Carray by distance
        struct dist_sort {
          inline bool operator()(Cpoint const& first, Cpoint const& second) const {
            return (first.dist_ < second.dist_);
          }
        };
        double Dist()         const { return dist_; }
        int PointsWithinEps() const { return pointsWithinEps_; }
        double Density() const { return density_; }
        int Fnum()       const { return fnum_; }
        int NearestIdx() const { return nidx_; }
        int Oidx()       const { return oidx_; }
        int Cnum()       const { return cnum_; }
        void SetPointsWithinEps(int d) { pointsWithinEps_ = d; }
        void SetDist(double d)         { dist_ = d; }
        void SetNearestIdx(int n)      { nidx_ = n; }
        void SetCluster(int c)         { cnum_ = c; }
        void AddDensity(double d)      { density_ += d; }
      private:
        double dist_; ///< minimum distance to point with higher density
        double density_; ///< Density from Gaussian kernel.
        int pointsWithinEps_; ///< # other points within epsilon
        int fnum_;    ///< Frame number.
        int nidx_;    ///< Index in Carray of nearest neighbor with higher density.
        int oidx_;    ///< Original index in Carray before sorting.
        int cnum_;    ///< Cluster number. -1 is no cluster.
    };
    // -------------------------------------------
    typedef std::vector<Cpoint> Carray;
    Carray Points_; ///< Hold info for each point to be clustered.
};

}
}
#endif

#ifndef INC_CLUSTER_DPEAKS_H
#define INC_CLUSTER_DPEAKS_H
#include "ClusterList.h"
class Cluster_DPeaks : public ClusterList {
  public:
    Cluster_DPeaks();
    static void Help();
    int SetupCluster(ArgList&);
    void ClusteringInfo();
    int Cluster();
    void AddSievedFrames();
    void ClusterResults(CpptrajFile&) const;
  private:
    void AssignClusterNum(int, int&);

    double epsilon_;
    class Cpoint {
      public:
        Cpoint() : dist_(-1.0), density_(0), fnum_(-1), nidx_(-1), cnum_(-1) {}
        Cpoint(int f) : dist_(-1.0), density_(0), fnum_(f), nidx_(-1), cnum_(-1) {}
        Cpoint(Cpoint const& rhs) : dist_(rhs.dist_), density_(rhs.density_), fnum_(rhs.fnum_),
                                    nidx_(rhs.nidx_), cnum_(rhs.cnum_) {}
        Cpoint& operator=(Cpoint const& rhs) {
          if (&rhs != this) {
            dist_ = rhs.dist_; density_ = rhs.density_; fnum_ = rhs.fnum_;
            nidx_ = rhs.nidx_; cnum_ = rhs.cnum_;
          }
          return *this;
        }
        /// Used to sort Carray
        bool operator<(Cpoint const& second) const {
          if (density_ == second.density_)
            return (fnum_ < second.fnum_);
          else
            return (density_ < second.density_);
        }
        double Dist()    const { return dist_; }
        int Density()    const { return density_; }
        int Fnum()       const { return fnum_; }
        int NearestIdx() const { return nidx_; }
        int Cnum()       const { return cnum_; }
        void SetDensity(int d)    { density_ = d; }
        void SetDist(double d)    { dist_ = d; }
        void SetNearestIdx(int n) { nidx_ = n; }
        void SetCluster(int c)    { cnum_ = c; }
      private:
        double dist_; ///< minimum distance to point with higher density
        int density_; ///< # other points within epsilon
        int fnum_;    ///< Frame number.
        int nidx_;    ///< Index in Carray of nearest neighbor with higher density.
        int cnum_;    ///< Cluster number. -1 is no cluster.
    };
    typedef std::vector<Cpoint> Carray;
    Carray Points_; ///< Hold info for each point to be clustered.
};
#endif

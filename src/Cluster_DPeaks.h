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
    double epsilon_;
    class Cpoint {
      public:
        Cpoint() : dist_(-1.0), density_(0), fnum_(-1) {}
        Cpoint(int f) : dist_(-1.0), density_(0), fnum_(f) {}
        /// Used to sort Carray
        bool operator<(Cpoint const& second) const {
          if (density_ == second.density_)
            return (fnum_ < second.fnum_);
          else
            return (density_ < second.density_);
        }
        double Dist() const { return dist_; }
        int Density() const { return density_; }
        int Fnum()    const { return fnum_; }
        void SetDensity(int d) { density_ = d; }
        void SetDist(double d) { dist_ = d; }
      private:
        double dist_; ///< minimum distance to point with higher density
        int density_; ///< # other points within epsilon
        int fnum_;    ///< Frame number.
    };
    typedef std::vector<Cpoint> Carray;
};
#endif

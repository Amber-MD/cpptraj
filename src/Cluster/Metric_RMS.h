#ifndef INC_CLUSTER_METRIC_RMS_H
#define INC_CLUSTER_METRIC_RMS_H
#include "Metric.h"
#include "../AtomMask.h"
#include "../DataSet_Coords.h"
#include "../Frame.h"
namespace Cpptraj {
namespace Cluster {

/// RMS cluster distance calc for Coords DataSet
class Metric_RMS : public Metric {
  public:
    Metric_RMS() : Metric(RMS), coords_(0), nofit_(false), useMass_(false) {}
    int Setup();
    double FrameDist(int, int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
    Metric* Copy() { return new Metric_RMS( *this ); }
    void FrameOpCentroid(int, Centroid*, double, CentOpType);
    std::string Description() const;
    void Info() const;
    unsigned int Ntotal() const { return (unsigned int)coords_->Size(); }
    // -------------------------------------------
    int Init(DataSet_Coords*, AtomMask const&, bool, bool);
    /// \return whether RMS is mass-weighted
    bool UseMass() const { return useMass_; }
    /// \return Atom mask
    AtomMask const& Mask() const { return mask_; }
  private:
    DataSet_Coords* coords_;
    AtomMask mask_;
    bool nofit_;
    bool useMass_;
    Frame frm1_;
    Frame frm2_;
};

}
}
#endif

#ifndef INC_CLUSTER_METRIC_QUATRMSD_H
#define INC_CLUSTER_METRIC_QUATRMSD_H
#include "Metric.h"
#include "../AtomMask.h"
#include "../DataSet_Coords.h"
#include "../Frame.h"
namespace Cpptraj {
namespace Cluster {

/// RMS cluster distance calc for Coords DataSet
class Metric_QuatRMSD : public Metric {
  public:
    Metric_QuatRMSD() : Metric(QRMSD), coords_(0), useMass_(false) {}
    int Setup();
    double FrameDist(int, int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
    Metric* Copy() { return new Metric_QuatRMSD( *this ); }
    void FrameOpCentroid(int, Centroid*, double, CentOpType);
    std::string Description() const;
    void Info() const;
    unsigned int Ntotal() const { return (unsigned int)coords_->Size(); }
    // -------------------------------------------
    int Init(DataSet_Coords*, AtomMask const&, bool);
    /// \return whether RMS is mass-weighted
    bool UseMass() const { return useMass_; }
    /// \return Atom mask
    AtomMask const& Mask() const { return mask_; }
  private:
    DataSet_Coords* coords_;
    AtomMask mask_;
    //bool nofit_;
    bool useMass_;
    Frame frm1_;
    Frame frm2_;
};

}
}
#endif

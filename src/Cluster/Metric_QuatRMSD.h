#ifndef INC_CLUSTER_METRIC_QUATRMSD_H
#define INC_CLUSTER_METRIC_QUATRMSD_H
#include "Metric.h"
#include "../AtomMask.h"
#include "../DataSet_Coords.h"
namespace Cpptraj {
namespace Cluster {

/// Symmetry-corrected RMS distance calc for Coords DataSet.
class Metric_QuatRMSD : public Metric {
  public:
    Metric_QuatRMSD() : Metric(QRMSD), useMass_(false) {}
    int Init(DataSet_Coords*,AtomMask const&,bool,bool,int);
    /// \return whether RMSD is mass-weighted
    bool UseMass() const { return useMass_; }
    /// \return Atom mask
    AtomMask const& Mask() const { return mask_; }
    // ----- Metric ------------------------------
    int Setup();
    double FrameDist(int, int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
    void FrameOpCentroid(int, Centroid*, double, CentOpType);
    Metric* Copy() { return new Metric_QuatRMSD( * this ); }
    std::string Description() const;
    void Info() const;
    unsigned int Ntotal() const { return (unsigned int)coords_->Size(); }
  private:
    DataSet_Coords* coords_;
    AtomMask mask_;
    Frame frm1_;
    Frame frm2_;
    bool useMass_;
}; 

}
}
#endif

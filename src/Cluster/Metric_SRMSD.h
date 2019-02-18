#ifndef INC_CLUSTER_METRIC_SRMSD_H
#define INC_CLUSTER_METRIC_SRMSD_H
#include "Metric.h"
#include "../AtomMask.h"
#include "../DataSet_Coords.h"
#include "../SymmetricRmsdCalc.h"
namespace Cpptraj {
namespace Cluster {

/// Symmetry-corrected RMS distance calc for Coords DataSet.
class Metric_SRMSD : public Metric {
  public:
    Metric_SRMSD() : Metric(SRMSD) {}
    int Init(DataSet_Coords*,AtomMask const&,bool,bool,int);
    // ----- Metric ------------------------------
    int Setup();
    double FrameDist(int, int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
    void FrameOpCentroid(int, Centroid*, double, CentOpType);
    Metric* Copy() { return new Metric_SRMSD( * this ); }
    std::string Description() const;
    void Info() const;
    unsigned int Ntotal() const { return (unsigned int)coords_->Size(); }
  private:
    DataSet_Coords* coords_;
    AtomMask mask_;
    SymmetricRmsdCalc SRMSD_;
    Frame frm1_;
    Frame frm2_;
}; 

}
}
#endif

#ifndef INC_CLUSTER_METRIC_DME_H
#define INC_CLUSTER_METRIC_DME_H
#include "../AtomMask.h"
#include "../DataSet_Coords.h"
#include "Metric.h"
namespace Cpptraj {
namespace Cluster {

/// DME cluster distance calc for Coords DataSet
class Metric_DME : public Metric {
  public:
    Metric_DME() : Metric(DME), coords_(0) {}
    int Setup();
    double FrameDist(int, int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
    Metric* Copy() { return new Metric_DME( *this ); }
    void FrameOpCentroid(int, Centroid*, double, CentOpType);
    std::string Description() const;
    void Info() const;
    unsigned int Ntotal() const { return (unsigned int)coords_->Size(); }
    // -------------------------------------------
    int Init(DataSet_Coords*, AtomMask const&);
  private:
    DataSet_Coords* coords_;
    AtomMask mask_;
    Frame frm1_;
    Frame frm2_;
};

}
}
#endif

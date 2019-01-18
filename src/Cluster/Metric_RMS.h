#ifndef INC_CLUSTER_METRIC_RMS_H
#define INC_CLUSTER_METRIC_RMS_H
#include "../AtomMask.h"
#include "../DataSet_Coords.h"
#include "Metric.h"
namespace Cpptraj {
namespace Cluster {

/// RMS cluster distance calc for Coords DataSet
class Metric_RMS : public Metric {
  public:
    Metric_RMS() : Metric(RMS), coords_(0), nofit_(false), useMass_(false) {}
    double FrameDist(int, int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
    Metric* Copy() { return new Metric_RMS( *this ); }
    void FrameOpCentroid(int, Centroid*, double, CentOpType);
    std::string Description() const;
    unsigned int Ntotal() const { return (unsigned int)coords_->Size(); }
    // -------------------------------------------
    int Setup(DataSet_Coords*, AtomMask const&, bool, bool);
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

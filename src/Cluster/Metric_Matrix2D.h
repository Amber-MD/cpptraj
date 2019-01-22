#ifndef INC_CLUSTER_METRIC_MATRIX2D_H
#define INC_CLUSTER_METRIC_MATRIX2D_H
#include "Metric.h"
#include "../DataSet_2D.h"
namespace Cpptraj {
namespace Cluster {

/// Metric for distance between two points in a matrix
class Metric_Matrix2D : public Metric {
  public:
    Metric_Matrix2D() : Metric(MATRIX2D), matrix_(0) {}
    int Setup();
    double FrameDist(int, int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
    Metric* Copy() { return new Metric_Matrix2D( *this ); }
    void FrameOpCentroid(int, Centroid*, double, CentOpType);
    std::string Description() const;
    void Info() const;
    unsigned int Ntotal() const { return (unsigned int)matrix_->Size(); }
    // -------------------------------------------
    int Init(DataSet_2D*);
  private:
    DataSet_2D* matrix_;
};

}
}
#endif

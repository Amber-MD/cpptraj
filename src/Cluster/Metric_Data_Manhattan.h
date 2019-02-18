#ifndef INC_CLUSTER_METRIC_DATA_MANHATTAN
#define INC_CLUSTER_METRIC_DATA_MANHATTAN
#include "Metric_Data.h"
namespace Cpptraj {
namespace Cluster {

/// Metric for scalar DataSet(s), Manhattan distance.
class Metric_Data_Manhattan : public Metric_Data {
  public:
    Metric_Data_Manhattan() : Metric_Data(MANHATTAN) {}
    // ----- Metric ------------------------------
    double FrameDist(int, int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    Centroid* NewCentroid(Cframes const&);
    Metric* Copy() { return new Metric_Data_Manhattan( *this ); }
    std::string Description() const;
    void Info() const;
};

}
}
#endif

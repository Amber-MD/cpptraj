#ifndef INC_CLUSTER_METRIC_DATA_EUCLID
#define INC_CLUSTER_METRIC_DATA_EUCLID
#include "Metric_Data.h"
namespace Cpptraj {
namespace Cluster {

/// Metric for scalar DataSet(s), Euclidean distance.
class Metric_Data_Euclid : public Metric_Data {
  public:
    Metric_Data_Euclid() : Metric_Data(EUCLID) {}
    // ----- Metric ------------------------------
    double FrameDist(int, int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    Metric* Copy() { return new Metric_Data_Euclid( *this ); }
    std::string Description() const;
    void Info() const;
};

}
}
#endif

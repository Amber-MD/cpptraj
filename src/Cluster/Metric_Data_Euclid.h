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
    int Setup();
    double FrameDist(int, int);
    double CentroidDist( Centroid*, Centroid* );
    double FrameCentroidDist(int, Centroid*);
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
    Metric* Copy() { return new Metric_Data_Euclid( *this ); }
    void FrameOpCentroid(int, Centroid*, double, CentOpType);
    std::string Description() const;
    void Info() const;
    unsigned int Ntotal() const;
    // -------------------------------------------
    int Init(DsArray const&);
  private:
    typedef std::vector<DataSet_1D*> D1Array;
    D1Array dsets_;
    typedef std::vector<DistCalc> DcArray;
    DcArray dcalcs_;
};

}
}
#endif

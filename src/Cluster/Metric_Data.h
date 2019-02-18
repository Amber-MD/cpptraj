#ifndef INC_CLUSTER_METRIC_DATA_H
#define INC_CLUSTER_METRIC_DATA_H
#include <vector>
#include "Metric.h"
#include "../DataSet_1D.h"
namespace Cpptraj {
namespace Cluster {

/// Abstract base class for Metrics that use scalar DataSets
class Metric_Data : public Metric {
  public:
    Metric_Data(Type t) : Metric(t), ntotal_(0) {}
    virtual ~Metric_Data() {}

    /// Input to Metric_Data
    typedef std::vector<DataSet*> DsArray;
    // ----- Metric ------------------------------
    int Setup();
    void CalculateCentroid(Centroid*, Cframes const&);
    Centroid* NewCentroid(Cframes const&);
    void FrameOpCentroid(int, Centroid*, double, CentOpType);
    unsigned int Ntotal() const { return ntotal_; }
    // -------------------------------------------
    int Init(DsArray const&);
  protected:
    typedef double (*DistCalc)(double,double);
    typedef std::vector<DataSet_1D*> D1Array;
    typedef std::vector<DistCalc> DcArray;

    static double DistCalc_Dih(double,double);
    static double DistCalc_Std(double,double);
    static double AvgCalc_Dih(DataSet_1D const&, Cframes const&, double&, double&);
    static double AvgCalc_Std(DataSet_1D const&, Cframes const&);
    static double DistCalc_FrameCentroid(double, double, bool, double, CentOpType, double&, double&);
    std::string SetNames(std::string const&) const;
    // TODO private
    D1Array dsets_;  ///< The input data sets
    DcArray dcalcs_; ///< Distance calc for each set
    unsigned int ntotal_; ///< Total number of data points in the smallest set
};

}
}
#endif

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
    Metric_Data(Type t) : Metric(t) {}
    virtual ~Metric_Data() {}

    /// Input to Metric_Data
    typedef std::vector<DataSet*> DsArray;
  protected:
    typedef double (*DistCalc)(double,double);

    static double DistCalc_Dih(double,double);
    static double DistCalc_Std(double,double);
    static double AvgCalc_Dih(DataSet_1D const&, Cframes const&, double&, double&);
    static double AvgCalc_Std(DataSet_1D const&, Cframes const&);
    static double DistCalc_FrameCentroid(double, double, bool, double, CentOpType, double&, double&);
};

}
}
#endif

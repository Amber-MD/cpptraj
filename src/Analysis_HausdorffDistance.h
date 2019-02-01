#ifndef INC_ANALYSIS_HAUSDORFFDISTANCE_H
#define INC_ANALYSIS_HAUSDORFFDISTANCE_H
#include "Analysis.h"
#include "DataSet_2D.h"
/// Compute the Hausdorff distance between 2 or more data sets.
/** The Hausdorff distance between two sets of points A and B is defined as:
  *   h(A,B) = max[a in A]{ min[b in B]( d(a,b) } }
  */
class Analysis_HausdorffDistance : public Analysis {
  public:
    Analysis_HausdorffDistance();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_HausdorffDistance(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    static double CalcHausdorffFromMatrix(DataSet_2D const&, double&, double&);

    enum OutType { BASIC = 0, UPPER_TRI_MATRIX, FULL_MATRIX };

    DataSetList inputSets_;
    OutType outType_;
    DataSet* out_;
    DataSet* ab_out_;
    DataSet* ba_out_;
};
#endif

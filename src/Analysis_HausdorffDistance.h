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
    Analysis_HausdorffDistance() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_HausdorffDistance(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    static double h_Matrix(DataSet_2D*);

    DataSetList inputSets_;
};
#endif

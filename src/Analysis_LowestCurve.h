#ifndef INC_ANALYSIS_LOWESTCURVE_H
#define INC_ANALYSIS_LOWESTCURVE_H
#include "Analysis.h"
#include "Array1D.h"
class Analysis_LowestCurve : public Analysis {
  public:
    Analysis_LowestCurve();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_LowestCurve(); }
    void Help() const;
  
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    Array1D input_dsets_;
    typedef std::vector<DataSet*> OutArray;
    OutArray output_sets_;
    int points_;
    double step_;
};
#endif

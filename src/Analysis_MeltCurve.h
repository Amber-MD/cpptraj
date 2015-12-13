#ifndef INC_ANALYSIS_MELTCURVE_H
#define INC_ANALYSIS_MELTCURVE_H
#include "Analysis.h"
#include "Array1D.h"
class Analysis_MeltCurve : public Analysis {
  public:
    Analysis_MeltCurve() : mcurve_(0), cut_(0.0) {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_MeltCurve(); }
    void Help() const;
  
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    Array1D input_dsets_;
    DataSet* mcurve_;
    double cut_;
};
#endif

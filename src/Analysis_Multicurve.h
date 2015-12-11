#ifndef INC_ANALYSIS_MULTICURVE_H
#define INC_ANALYSIS_MULTICURVE_H
#include "Analysis.h"
#include "Array1D.h"
class Analysis_Multicurve : public Analysis {
  public:
    Analysis_Multicurve() : debug_(0) {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Multicurve(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    Array1D inputDsets_;
    ArgList args_;
    AnalysisSetup master_;
    int debug_;
};
#endif

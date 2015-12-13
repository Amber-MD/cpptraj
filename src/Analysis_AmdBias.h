#ifndef INC_ANALYSIS_AMDBIAS_H
#define INC_ANALYSIS_AMDBIAS_H
#include "Analysis.h"
class Analysis_AmdBias : public Analysis {
  public:
    Analysis_AmdBias();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_AmdBias(); }
    void Help() const;
  
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    DataSet* ds1_;
    double Ethresh_;
    double alpha_;
    DataSet* bias_;
};
#endif

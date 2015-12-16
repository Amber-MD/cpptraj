#ifndef INC_ANALYSIS_RUNNNINGAVG_H
#define INC_ANALYSIS_RUNNNINGAVG_H
#include "Analysis.h"
#include "Array1D.h"
class Analysis_RunningAvg : public Analysis {
  public:
    Analysis_RunningAvg();

    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_RunningAvg(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    Array1D dsets_;
    bool cumulative_;
    int window_;
    std::vector<DataSet*> outputData_;
};
#endif

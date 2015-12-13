#ifndef INC_ANALYSIS_AUTOCORR_H
#define INC_ANALYSIS_AUTOCORR_H
#include "Analysis.h"
class Analysis_AutoCorr : public Analysis {
  public:
    Analysis_AutoCorr();

    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_AutoCorr(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    DataSetList::DataListType dsets_;
    DataSetList::DataListType outputData_;
    int lagmax_;
    bool usefft_;
    bool calc_covar_;
};
#endif

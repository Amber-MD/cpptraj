#ifndef INC_ANALYSIS_AUTOCORR_H
#define INC_ANALYSIS_AUTOCORR_H
#include "Analysis.h"
class Analysis_AutoCorr : public Analysis {
  public:
    Analysis_AutoCorr();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_AutoCorr(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    DataSetList::DataListType dsets_;
    DataSetList::DataListType outputData_;
    int lagmax_;
    bool usefft_;
    bool calc_covar_;
};
#endif

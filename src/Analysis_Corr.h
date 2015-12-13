#ifndef INC_ANALYSIS_CORR_H
#define INC_ANALYSIS_CORR_H
#include "Analysis.h"
// Class: Analysis_Corr
/// Calculate autocorrelation or correlation of dataset(s).
class Analysis_Corr : public Analysis {
  public:
    Analysis_Corr();

    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Corr(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    DataSet *D1_;
    DataSet *D2_;
    int lagmax_;
    DataSet* Ct_;
    bool usefft_;
    bool calc_covar_;
};
#endif

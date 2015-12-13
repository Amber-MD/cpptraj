#ifndef INC_ANALYSIS_CROSSCORR_H
#define INC_ANALYSIS_CROSSCORR_H
#include "Analysis.h"
#include "Array1D.h"
class Analysis_CrossCorr : public Analysis {
  public:
    Analysis_CrossCorr();

    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_CrossCorr(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();

  private:
    DataFile* outfile_;
    Array1D input_dsets_;
    DataSet* matrix_;
};
#endif

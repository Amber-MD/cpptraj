#ifndef INC_ANALYSIS_MULTIHIST_H
#define INC_ANALYSIS_MULTIHIST_H
#include "Analysis.h"
/// Histogram multiple 1D data sets separately.
class Analysis_MultiHist : public Analysis {
  public:
    Analysis_MultiHist();
    ~Analysis_MultiHist();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_MultiHist(); }
    void Help() const;
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    typedef std::vector<Analysis*> Harray; // Use ptrs, no copy construct defined
    Harray Histograms_;
};
#endif

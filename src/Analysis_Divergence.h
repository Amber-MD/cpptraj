#ifndef INC_ANALYSIS_DIVERGENCE_H
#define INC_ANALYSIS_DIVERGENCE_H
#include "Analysis.h"
class Analysis_Divergence : public Analysis {
  public:
    Analysis_Divergence();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Divergence(); }
    void Help() const;
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    /// Normalize sum over set to 1.0
    std::vector<double> NormalizeSet(DataSet const&, unsigned int) const;
    DataSet* ds1_;
    DataSet* ds2_;
    DataSet* dataout_;
};
#endif

#ifndef INC_ANALYSIS_CONSTANTPHSTATS_H
#define INC_ANALYSIS_CONSTANTPHSTATS_H
#include "Analysis.h"
/// <Enter description of Analysis_ConstantPHStats here>
class Analysis_ConstantPHStats : public Analysis {
  public:
    Analysis_ConstantPHStats() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_ConstantPHStats(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    int debug_;
    DataSetList inputSets_;
};
#endif

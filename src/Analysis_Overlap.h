#ifndef INC_ANALYSIS_OVERLAP_H
#define INC_ANALYSIS_OVERLAP_H
#include "Analysis.h"
class Analysis_Overlap : public Analysis {
  public:
    Analysis_Overlap();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Overlap(); }
    void Help() const;
  
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    DataSet* ds1_;
    DataSet* ds2_;
    bool useDeviation_;
};
#endif

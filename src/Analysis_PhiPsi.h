#ifndef INC_ANALYSIS_PHIPSI_H
#define INC_ANALYSIS_PHIPSI_H
#include "Analysis.h"
#include "Array1D.h"
class Analysis_PhiPsi : public Analysis {
  public:
    Analysis_PhiPsi();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_PhiPsi(); }
    void Help() const;
  
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    Array1D input_dsets_;
    CpptrajFile outfile_;
};
#endif

#ifndef INC_ANALYSIS_EVALEQUILIBRATION_H
#define INC_ANALYSIS_EVALEQUILIBRATION_H
#include "Analysis.h"
#include "Array1D.h"
/// <Enter description of Analysis_EvalEquilibration here>
class Analysis_EvalEquilibration : public Analysis {
  public:
    Analysis_EvalEquilibration();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_EvalEquilibration(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    Array1D inputSets_;
    std::string dsname_; ///< Output set name
    int debug_;
};
#endif

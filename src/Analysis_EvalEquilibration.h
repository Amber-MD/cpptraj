#ifndef INC_ANALYSIS_EVALEQUILIBRATION_H
#define INC_ANALYSIS_EVALEQUILIBRATION_H
#include "Analysis.h"
/// <Enter description of Analysis_EvalEquilibration here>
class Analysis_EvalEquilibration : public Analysis {
  public:
    Analysis_EvalEquilibration();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_EvalEquilibration(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    DataSet* setIn_;
    std::string dsname_; ///< Output set name
    int debug_;
};
#endif

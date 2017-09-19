#ifndef INC_ANALYSIS_REGRESSION_H
#define INC_ANALYSIS_REGRESSION_H
#include "Analysis.h"
#include "Array1D.h"
class Analysis_Regression : public Analysis {
  public:
    Analysis_Regression();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Regression(); }
    void Help() const;
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    Array1D input_dsets_;
    Array1D output_dsets_;
    typedef std::vector<DataSet*> DSarray;
    DSarray slope_dsets_;
    DSarray int_dsets_;
    CpptrajFile* statsout_;
    int nx_; ///< Number of x values
};
#endif

#ifndef INC_ANALYSIS_INTEGRATE_H
#define INC_ANALYSIS_INTEGRATE_H
#include "Analysis.h"
#include "Array1D.h"
class DataSet_Mesh;
class Analysis_Integrate : public Analysis {
  public:
    Analysis_Integrate();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Integrate(); }
    void Help() const;
  
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    Array1D input_dsets_;
    std::vector<DataSet_Mesh*> output_dsets_;
    DataSet* sumSet_; ///< Hold final sum for each set
};
#endif

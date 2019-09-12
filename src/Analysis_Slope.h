#ifndef INC_ANALYSIS_SLOPE_H
#define INC_ANALYSIS_SLOPE_H
#include "Analysis.h"
#include "Array1D.h"
class DataSet_Mesh;
/// Calculate the slope (finite difference) for input DataSets
class Analysis_Slope : public Analysis {
  public:
    Analysis_Slope() : diffType_(DataSet_1D::FORWARD) {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Slope(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    DataSet_1D::DiffType diffType_;
    Array1D input_dsets_;
    std::vector<DataSet_Mesh*> output_dsets_;
};
#endif

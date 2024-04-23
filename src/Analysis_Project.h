#ifndef INC_ANALYSIS_PROJECT_H
#define INC_ANALYSIS_PROJECT_H
#include "Analysis.h"
#include "Array1D.h"
class DataSet_Modes;
/// <Enter description of Analysis_Project here>
class Analysis_Project : public Analysis {
  public:
    Analysis_Project() {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Project(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    int beg_;
    int end_;
    DataSet_Modes* modinfo_; ///< Modes to project on
    Array1D Sets_; ///< Input data sets to project
    std::vector<DataSet*> project_; ///< Hold output projections for each mode
};
#endif

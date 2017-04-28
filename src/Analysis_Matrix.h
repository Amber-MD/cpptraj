#ifndef INC_ANALYSIS_MATRIX_H
#define INC_ANALYSIS_MATRIX_H
#include "Analysis.h"
#include "DataSet_2D.h"
#include "DataSet_Modes.h"
class Analysis_Matrix : public Analysis {
  public:
    Analysis_Matrix();

    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Matrix(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    int NMWizOutput() const;

    DataSet_2D* matrix_;
    DataSet_Modes* modes_;
    CpptrajFile* outthermo_;
    double thermo_temp_;
    int nevec_;
    bool thermopt_;
    bool reduce_;
    bool nmwizopt_;
    int nmwizvecs_;
    CpptrajFile* nmwizfile_;
    Topology nmwizParm_;
};
#endif

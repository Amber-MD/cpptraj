#ifndef INC_ANALYSIS_CLUSTER_H
#define INC_ANALYSIS_CLUSTER_H
#include "Analysis.h"
#include "Cluster/Control.h"
#include "DataSet_Coords.h"
/// <Enter description of Analysis_Cluster here>
class Analysis_Cluster : public Analysis {
  public:
    Analysis_Cluster() : coords_(0), masterDSL_(0) {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Cluster(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    Cpptraj::Cluster::Control control_;

    DataSet_Coords* coords_;
    DataSetList* masterDSL_;
};
#endif

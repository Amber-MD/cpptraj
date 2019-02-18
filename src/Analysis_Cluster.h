#ifndef INC_ANALYSIS_CLUSTER_H
#define INC_ANALYSIS_CLUSTER_H
#include "Analysis.h"
#include "Cluster/Control.h"
/// <Enter description of Analysis_Cluster here>
class Analysis_Cluster : public Analysis {
  public:
    Analysis_Cluster() : masterDSL_(0) {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Cluster(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    Cpptraj::Cluster::Control control_;

    DataSetList* masterDSL_;
};
#endif

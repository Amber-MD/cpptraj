#ifndef INC_ANALYSIS_CLUSTERING_H
#define INC_ANALYSIS_CLUSTERING_H
#include "Analysis.h"
#include "Cluster/Control.h"
/// Perform cluster analysis on DataSets 
class Analysis_Clustering : public Analysis {
  public:
    Analysis_Clustering() : masterDSL_(0) {}
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_Clustering(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    Cpptraj::Cluster::Control control_;

    DataSetList* masterDSL_;
};
#endif

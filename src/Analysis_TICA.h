#ifndef INC_ANALYSIS_TICA_H
#define INC_ANALYSIS_TICA_H
#include "Analysis.h"
/// <Enter description of Analysis_TICA here>
class Analysis_TICA : public Analysis {
  public:
    Analysis_TICA();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_TICA(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    DataSet_Coords* TgtTraj_; ///< Input trajectory
    int lag_; ///< TICA time lag
    AtomMask mask1_; ///< Atoms to use in matrix calc
};
#endif

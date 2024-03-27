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
    AtomMask mask2_; ///< Second atom mask for debugging full covar matrix
    bool useMass_; ///< Control whether to mass-weight
    CpptrajFile* debugFile_; ///< Debug output
};
#endif

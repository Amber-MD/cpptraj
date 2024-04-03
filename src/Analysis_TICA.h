#ifndef INC_ANALYSIS_TICA_H
#define INC_ANALYSIS_TICA_H
#include "Analysis.h"
#include "Array1D.h"
/// <Enter description of Analysis_TICA here>
class Analysis_TICA : public Analysis {
  public:
    Analysis_TICA();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_TICA(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    /// Analyze using coordinates data set (TgtTraj_)
    Analysis::RetType analyze_crdset();
    /// Analyze using 1D data sets (sets_)
    Analysis::RetType analyze_datasets();

    DataSet_Coords* TgtTraj_; ///< Input trajectory (crdset)
    int lag_; ///< TICA time lag
    AtomMask mask1_; ///< Atoms to use in matrix calc
    AtomMask mask2_; ///< Second atom mask for debugging full covar matrix
    bool useMass_; ///< Control whether to mass-weight
    CpptrajFile* debugC0_; ///< Debug output for C0
    CpptrajFile* debugCT_; ///< Debug output for CT
    Array1D sets_;         ///< 1D data sets (data)
};
#endif

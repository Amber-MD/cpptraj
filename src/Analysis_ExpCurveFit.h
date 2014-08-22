#ifndef INC_ANALYSIS_EXPCURVEFIT_H
#define INC_ANALYSIS_EXPCURVEFIT_H
#include "Analysis.h"
/// Perform curve fit to sum of exponentials.
class Analysis_ExpCurveFit : public Analysis {
  public:
    Analysis_ExpCurveFit();
                          
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_ExpCurveFit(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    enum EqFormType { MEXP, MEXP_K, MEXP_K_PENALTY }; 

    DataSet* dset_;
    DataSet* finalY_;
    double tolerance_;
    int maxIt_;
    int nexp_;          ///< Number of exponentials
    EqFormType eqForm_; ///< Equation form
    bool usePrefactorBounds_; 
};
#endif

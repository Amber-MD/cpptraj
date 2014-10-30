#ifndef INC_ANALYSIS_CURVEFIT
#define INC_ANALYSIS_CURVEFIT
#include "Analysis.h"
/// Used for non-linear curve fitting.
class Analysis_CurveFit : public Analysis {
  public:
    Analysis_CurveFit();
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_CurveFit(); }
    static void Help();
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    enum EqFormType { GENERAL = 0, MEXP, MEXP_K, MEXP_K_PENALTY };
    DataSet* dset_;     ///< DataSet to fit.
    DataSet* finalY_;   ///< Final output DataSet.
    typedef std::vector<double> Darray;
    Darray Params_;     ///< Equation parameters.
    double tolerance_;  ///< Curve fit tolerance.
    int maxIt_;         ///< Max # iterations.
    int nexp_;          ///< # exponentials.
    EqFormType eqForm_; ///< Equation form.
};
#endif

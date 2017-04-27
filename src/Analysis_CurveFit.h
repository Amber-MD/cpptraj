#ifndef INC_ANALYSIS_CURVEFIT
#define INC_ANALYSIS_CURVEFIT
#include "Analysis.h"
/// Used for non-linear curve fitting.
class Analysis_CurveFit : public Analysis {
  public:
    Analysis_CurveFit();
    Analysis_CurveFit(DataSet*, int, ArgList&, DataSetList&, DataFileList&, int);

    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_CurveFit(); }
    static void HelpText(); ///< For use in Analysis_MultiCurve
    void Help() const;
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    Analysis::RetType Internal_setup(std::string const&, ArgList&, DataSetList&, DataFileList&, int);

    enum EqFormType { GENERAL = 0, MEXP, MEXP_K, MEXP_K_PENALTY, GAUSS };
    std::string equation_; ///< Equation to fit.
    CpptrajFile* Results_; ///< Results output file (final params, stats)
    DataSet* dset_;     ///< DataSet to fit.
    DataSet* finalY_;   ///< Final output DataSet.
    typedef std::vector<double> Darray;
    std::vector<DataSet*> A_param_sets_; ///< Hold final equation parameters.
    DataSet* set_corr_; ///< Hold correlation coefficient.
    DataSet* set_chi_;  ///< Hold chi-squared.
    DataSet* set_unc_;  ///< Hold uncertainty.
    DataSet* set_rms_;  ///< Hold RMS percent error.
    Darray Params_;     ///< Equation parameters.
    double tolerance_;  ///< Curve fit tolerance.
    double outXmin_;    ///< Output X min.
    double outXmax_;    ///< Output X max.
    int maxIt_;         ///< Max # iterations.
    int nexp_;          ///< # exponentials.
    int outXbins_;      ///< # of points in output DataSet.
    int n_expected_params_; ///< Number of expected parameters.
    int n_specified_params_;///< Number of specified parameters.
    EqFormType eqForm_; ///< Equation form.
};
#endif

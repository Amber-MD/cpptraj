#ifndef INC_ANALYSIS_KDE_H
#define INC_ANALYSIS_KDE_H
#include "Analysis.h"
#include "DataSet_1D.h"
class Analysis_KDE : public Analysis {
  public:
    Analysis_KDE();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_KDE(); }
    void Help() const;

    Analysis::RetType ExternalSetup(DataSet_1D*, std::string const&, int, std::string const&,
                            bool, double, bool, double, double, int, double,
                            DataSetList&, DataFileList&);
    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    static const double ONE_OVER_ROOT_TWOPI;
    typedef double (Analysis_KDE::*fxnptr)(double) const;

    double GaussianKernel(double) const;

    DataSet* data_;    ///< Data set to histogram.
    DataSet* q_data_;  ///< Second set if calculating KL divergence.
    double bandwidth_; ///< Bandwidth for KDE.
    DataSet* output_;  ///< Output Histogram.
    DataSet* kldiv_;   ///< KL divergence vs time.
    DataSet* amddata_; ///< Optional AMD boost data set.
    bool calcFreeE_;
    double Temp_;
    fxnptr Kernel_;    ///< Kernel to use.
    double default_min_;
    double default_max_;
    double default_step_;
    int default_bins_;
    bool minArgSet_;
    bool maxArgSet_;
};
#endif

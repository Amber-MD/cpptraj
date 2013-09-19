#ifndef INC_ANALYSIS_KDE_H
#define INC_ANALYSIS_KDE_H
#include "Analysis.h"
class Analysis_KDE : public Analysis {
  public:
    Analysis_KDE();
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_KDE(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
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
    fxnptr Kernel_;    ///< Kernel to use.
};
#endif

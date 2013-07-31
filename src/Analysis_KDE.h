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

    DataSet* data_;
    double bandwidth_;
    DataSet* output_;
    fxnptr Kernel_;
};
#endif

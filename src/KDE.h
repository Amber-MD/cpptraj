#ifndef INC_KDE_H
#define INC_KDE_H
#include "HistBin.h"
#include "DataSet_double.h"
/// Can be used to calculate 1D histogram using kernel density estimator.
class KDE {
  public:
    enum KernelType { GAUSSIAN = 0 };
    KDE();
#   ifdef _OPENMP
    /// CONSTRUCTOR - Number of threads.
    KDE(int);
#   endif
    /// Output, Input
    int CalcKDE(DataSet_double&, DataSet_1D const&) const;
    /// Output, Input, Increments, Histogram dimension, Bandwidth
    int CalcKDE(DataSet_double&, DataSet_1D const&,
                std::vector<double> const&,HistBin const&,double) const;
  private:
    static const double ONE_OVER_ROOT_TWOPI;
    typedef double (KDE::*FxnPtr)(double) const;
    double GaussianKernel(double) const;
#   ifdef _OPENMP
    int numthreads_;
#   endif
    KernelType ktype_;
    FxnPtr Kernel_;
};
#endif

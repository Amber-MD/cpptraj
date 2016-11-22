#ifndef INC_KDE_H
#define INC_KDE_H
#include "HistBin.h"
#include "DataSet_double.h"
/// Can be used to calculate 1D histogram using kernel density estimator.
class KDE {
  public:
    enum KernelType { NONE = 0, GAUSSIAN };
    KDE();
    int CalcKDE(DataSet_double&, DataSet_1D const&) const;
  private:
    static const double ONE_OVER_ROOT_TWOPI;
    typedef double (KDE::*FxnPtr)(double) const;
    double GaussianKernel(double) const;

    KernelType ktype_;
    FxnPtr Kernel_;
    
    /// Output, Input, Histogram dimension, Bandwidth
    int CalcKDE(DataSet_double&, DataSet_1D const&,
                std::vector<double> const&,HistBin const&,double) const;
};
#endif

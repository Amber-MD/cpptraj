#ifndef INC_KDE_H
#define INC_KDE_H
#include "HistBin.h"
#include "DataSet_double.h"
/// Can be used to calculate 1D histogram using kernel density estimator.
namespace KDE {
  typedef double (*FxnPtr)(double);
  /// Output, Input, Histogram dimension, Bandwidth
  int GaussianKDE(DataSet_double&,DataSet_1D const&,
                  std::vector<double> const&,HistBin const&,double);
  /// Kernel, Output, Input, Histogram dimension, Bandwidth
  int CalcKDE(FxnPtr, DataSet_double&,DataSet_1D const&,
              std::vector<double> const&,HistBin const&,double);
}
#endif

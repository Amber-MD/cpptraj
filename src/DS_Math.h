#ifndef INC_DS_MATH_H
#define INC_DS_MATH_H
/*! \file DS_Math.h
    \brief Collection of routines to perform math on 1D datasets.
 */
#include "DataSet_1D.h"
namespace DS_Math {
  double Avg(DataSet_1D const&, double*);
  double Avg(DataSet_1D const&);
  double Min(DataSet_1D const&);
  double Max(DataSet_1D const&);
  int CrossCorr(DataSet_1D const&, DataSet_1D const&, DataSet_1D&, int, bool, bool);
  double CorrCoeff(DataSet_1D const&, DataSet_1D const&);
}
#endif

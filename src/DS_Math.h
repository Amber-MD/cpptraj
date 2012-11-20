#ifndef INC_DS_MATH_H
#define INC_DS_MATH_H
/*! \file DS_Math.h
    \brief Collection of routines to perform math on datasets.
 */
#include "DataSet.h"
namespace DS_Math {
  // TODO: Make const refs
  double Avg(DataSet&, double*);
  double Avg(DataSet&);
  double Min(DataSet&);
  double Max(DataSet&);
  int CrossCorr(DataSet&, DataSet&, DataSet&, int, bool, bool);
  double CorrCoeff(DataSet&, DataSet&);
}
#endif

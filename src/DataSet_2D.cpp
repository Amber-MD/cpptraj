#include "DataSet_2D.h"
#include <cmath> // fabs

/** Check that the given matrix is symmetric. 
  * This will work similar to numpy.allclose.
  * https://numpy.org/doc/stable/reference/generated/numpy.allclose.html
  */
bool DataSet_2D::IsSymmetric() const {
  if (MatrixKind() == HALF || MatrixKind() == TRI) return true;
  if (Ncols() != Nrows()) return false;

  static const double rtol = 1e-5;
  static const double atol = 1e-8;

  for (unsigned int row = 0; row != Nrows(); row++) {
    for (unsigned int col = row+1; col < Ncols(); col++) {
      double a = GetElement(col, row);
      double b = GetElement(row, col);
      bool is_equiv = ((fabs(a - b) <= (atol + rtol * fabs(b))));
      if (!is_equiv) return false;
    }
  }
  return true;
}

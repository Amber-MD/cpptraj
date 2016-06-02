#include <cfloat> // FLT_MAX
#include "ClusterMatrix.h"

// ClusterMatrix::FindMin()
/** Find the minimum; set corresponding row and column. Cannot currently
  * be used for sieved frames.
  */
double ClusterMatrix::FindMin(int& iOut, int& jOut) const {
  float min = FLT_MAX;
  for (unsigned int row = 0; row != Mat_.Nrows(); row++) {
    if (!ignore_[row]) {
      unsigned int col = row + 1;
      unsigned int idx = Mat_.CalcIndex(col, row); // idx is start of this row
      for (; col != Mat_.Ncols(); col++, idx++) {
        if (!ignore_[col] && Mat_[idx] < min) {
          min = Mat_[idx];
          iOut = (int)row;
          jOut = (int)col;
        }
      }
    }
  }
  return (double)min;
}

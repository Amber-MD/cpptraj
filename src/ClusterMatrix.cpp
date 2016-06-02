#include <cfloat> // FLT_MAX
#include "ClusterMatrix.h"
#include "CpptrajStdio.h" // PrintElements

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

void ClusterMatrix::PrintElements() const {
  unsigned int iVal = 0;
  unsigned int jVal = 1;
  for (size_t idx = 0UL; idx < Mat_.size(); ++idx) {
    if (!ignore_[iVal] && !ignore_[jVal])
      mprintf("\t%u %u %f\n",iVal,jVal,Mat_[idx]);
    // Increment indices
    jVal++;
    if (jVal >= ignore_.size()) {
      iVal++;
      jVal = iVal + 1;
    }
  }
}

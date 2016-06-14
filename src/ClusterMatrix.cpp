#include <cfloat> // FLT_MAX
#include "ClusterMatrix.h"
#include "CpptrajStdio.h" // PrintElements
#ifdef _OPENMP
# include <omp.h>
#endif

// ClusterMatrix::FindMin()
/** Find the minimum; set corresponding row and column. Cannot currently
  * be used for sieved frames.
  */
double ClusterMatrix::FindMin(int& iOut, int& jOut) const {
# ifdef _OPENMP
  static int minRow_[128]; // FIXME should be allocd eventually
  static int minCol_[128];
  static float minVal_[128];
  int row, mythread, numthreads;
  int nrows = (int)Mat_.Nrows();
# pragma omp parallel private(row, mythread)
  {
  mythread = omp_get_thread_num();
  if (mythread == 0) numthreads = omp_get_num_threads();
  minVal_[mythread] = FLT_MAX;
# pragma omp for schedule(dynamic)
  for (row = 0; row < nrows; row++) {
    if (!ignore_[row]) {
      int col = row + 1;
      unsigned int idx = Mat_.CalcIndex(col, row); /// idx is start of this row
      for (; col != nrows; col++, idx++) {
        if (!ignore_[col]) {
          if (Mat_[idx] < minVal_[mythread]) {
            minVal_[mythread] = Mat_[idx];
            minRow_[mythread] = row;
            minCol_[mythread] = col;
          }
        }
      }
    }
  }
  } /* END pragma omp parallel */
  float min = minVal_[0];
  iOut = minRow_[0];
  jOut = minCol_[0];
  for (int idx = 1; idx != numthreads; idx++) {
    if (minVal_[idx] < min) {
      min = minVal_[idx];
      iOut = minRow_[idx];
      jOut = minCol_[idx];
    }
  }
# else
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
# endif
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

// ClusterMatrix::SetupMatrix()
int ClusterMatrix::SetupMatrix(size_t sizeIn) {
  if (Mat_.resize( 0L, sizeIn )) return 1;
  ignore_.assign( sizeIn, false );
/*
# ifdef _OPENMP
  int n_threads = 0;
# pragma omp parallel
  {
    if (omp_get_thread_num() == 0)
      n_threads = omp_get_num_threads();
  }
  minRow_.resize( n_threads );
  minCol_.resize( n_threads );
  minVal_.resize( n_threads );
  mprintf("DEBUG: ClusterMatrix: Using %i threads.\n", n_threads);
# endif
*/
  return 0;
}

#include <limits> // float max
#ifdef _OPENMP
#include <omp.h>
#endif
#include "DynamicMatrix.h"
#include "../CpptrajStdio.h"

// DynamicMatrix::FindMin()
/** Find the minimum not being ignored; set corresponding row and column. */
#ifdef _OPENMP
double Cpptraj::Cluster::DynamicMatrix::FindMin(int& iOut, int& jOut) {
  int row, mythread;
  int nrows = (int)Mat_.Nrows();
# pragma omp parallel private(row, mythread)
  {
  mythread = omp_get_thread_num();

  minVal_[mythread] = std::numeric_limits<float>::max();
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
  for (unsigned int idx = 1; idx != minVal_.size(); idx++) {
    if (minVal_[idx] < min) {
      min = minVal_[idx];
      iOut = minRow_[idx];
      jOut = minCol_[idx];
    }
  }
  return (double)min;
}
#else
double Cpptraj::Cluster::DynamicMatrix::FindMin(int& iOut, int& jOut) const {
  float min = std::numeric_limits<float>::max();
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
#endif

void Cpptraj::Cluster::DynamicMatrix::PrintElements() const {
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

// DynamicMatrix::SetupMatrix()
int Cpptraj::Cluster::DynamicMatrix::SetupMatrix(size_t sizeIn) {
  if (Mat_.resize( 0L, sizeIn )) return 1;
  ignore_.assign( sizeIn, false );
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
# endif
  return 0;
}

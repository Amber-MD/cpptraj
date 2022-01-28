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
/*
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
  mprintf("DEBUG: Min found at row %i col %i (%f)\n", iOut, jOut, min);*/


  float currentMin = std::numeric_limits<float>::max();
  int minRow = -1;
  int minCol = -1;
  for (unsigned int col = 0; col != closestIdx_.size(); col++)
  {
    int row = closestIdx_[col];
    if (!ignore_[col] && !ignore_[row]) {
      float mval = Mat_.element(col, row);
      if (mval < currentMin) {
        currentMin = mval;
        minRow = (int)row;
        minCol = (int)col;
      }
    }
  }
  // Ensure iOut < jOut
  if (minRow < minCol) {
    iOut = minRow;
    jOut = minCol;
  } else {
    iOut = minCol;
    jOut = minRow;
  }
  return (double)currentMin;
/*
  if (minCol < minRow) { //FIXME i think only needed for testing
    int itmp = minCol;
    minCol = minRow;
    minRow = itmp;
  }
  mprintf("DEBUG: Second try: Min found at row %i col %i (%f)\n", minRow, minCol, currentMin);

  if (iOut != minRow || jOut != minCol)
    mprintf("Error: Column/row mismatch.\n");

  return (double)min;*/
}
#endif

/** For DEBUG. Print all non-ignored matrix elements. */
void Cpptraj::Cluster::DynamicMatrix::PrintElements() const {
  unsigned int iVal = 0;
  unsigned int jVal = 1;
  for (size_t idx = 0UL; idx < Mat_.size(); ++idx) {
    if (!ignore_[iVal] && !ignore_[jVal])
      mprintf("\t%u %u %20.10E\n",iVal,jVal,Mat_[idx]);
    // Increment indices
    jVal++;
    if (jVal >= ignore_.size()) {
      iVal++;
      jVal = iVal + 1;
    }
  }
}

/** For DEBUG. Print elements like a square. */
void Cpptraj::Cluster::DynamicMatrix::PrintElementsSquare() const {
/*  mprintf("%4s", "");
  for (unsigned int col = 0; col < Mat_.Ncols(); col++)
    mprintf(" %12u", col);
  mprintf("\n");
  for (unsigned int row = 0; row < Mat_.Nrows(); row++) {
    mprintf("%4u", row);
    for (unsigned int col = 0; col < Mat_.Ncols(); col++) {
      if (row <= col)
        mprintf(" %12s", "X");
      else if (ignore_[row] || ignore_[col])
        mprintf(" %12s", "i");
      else
        mprintf(" %12.4f", Mat_.element(col, row));
    }
    mprintf("\n");
  }*/
  for (unsigned int row = 0; row < Mat_.Nrows(); row++) {
    if (!ignore_[row]) {
      mprintf("%4u :", row);
      for (unsigned int col = 0; col < Mat_.Ncols(); col++) {
        if (col != row && !ignore_[col])
          mprintf(" %4u=%8.3f", col, Mat_.element(col, row));
      }
      mprintf("\n");
    }
  }
}

/** For DEBUG. Print closest indices to STDOUT. */
void Cpptraj::Cluster::DynamicMatrix::PrintClosest() const {
  for (unsigned int idx = 0; idx != closestIdx_.size(); idx++) {
    if (!ignore_[idx]) {
      mprintf("DEBUG: Closest cluster to %u is %i (%f)\n", idx, closestIdx_[idx], Mat_.element(idx, closestIdx_[idx]));
      // Check that this is actually the case
      int cidx = -1;
      float cdist = 0;
      for (unsigned int jdx = 0; jdx != closestIdx_.size(); jdx++) {
        if (idx != jdx && !ignore_[jdx]) {
          float fdist = Mat_.element(idx, jdx);
          if (cidx == -1) {
            cidx = jdx;
            cdist = fdist;
          } else if (fdist < cdist) {
            cidx = jdx;
            cdist = fdist;
          }
        }
      }
      if (cidx != closestIdx_[idx]) {
        mprintf("Error: ACTUAL closest idx to %u is %i (%f)\n", idx, cidx, cdist);
      }
    }
  }
}

// DynamicMatrix::SetupMatrix()
int Cpptraj::Cluster::DynamicMatrix::SetupMatrix(size_t sizeIn) {
  if (Mat_.resize( 0L, sizeIn )) return 1;
  ignore_.assign( sizeIn, false );
  closestIdx_.assign( sizeIn, -1 );
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

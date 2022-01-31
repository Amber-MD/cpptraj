#include <limits> // float max
#include "DynamicMatrix.h"
#include "../CpptrajStdio.h"

// DynamicMatrix::FindMin()
/** Find the minimum not being ignored; set corresponding row and column. */
double Cpptraj::Cluster::DynamicMatrix::FindMin(int& iOut, int& jOut) const {
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
}

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
  return 0;
}

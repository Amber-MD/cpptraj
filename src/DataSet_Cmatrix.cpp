#include <cfloat> // FLT_MAX
#include "DataSet_Cmatrix.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

void DataSet_Cmatrix::WriteBuffer(CpptrajFile& outfile, SizeArray const& pIn) const {
  size_t x = (size_t)pIn[0];
  size_t y = (size_t)pIn[1];
  if ( x >= Mat_.Ncols() || y >= Mat_.Nrows() )
    outfile.Printf(format_.fmt(), 0.0);
  else
    outfile.Printf(format_.fmt(), Mat_.element(x,y));
}

double* DataSet_Cmatrix::MatrixArray() const {
  double* matOut = new double[ Mat_.size() ];
  for (size_t i = 0; i < Mat_.size(); ++i)
    matOut[i] = (double)Mat_[i];
  return matOut;
}

// -----------------------------------------------------------------------------
// COPY CONSTRUCTOR
DataSet_Cmatrix::DataSet_Cmatrix(const DataSet_Cmatrix& rhs) :
  DataSet_2D( rhs ),
  ignore_(rhs.ignore_),
  Mat_(rhs.Mat_)
{}

// ASSIGNMENT
DataSet_Cmatrix& DataSet_Cmatrix::operator=(const DataSet_Cmatrix& rhs) {
  if (this != &rhs) {
    DataSet_2D::operator=( rhs );
    ignore_ = rhs.ignore_;
    Mat_ = rhs.Mat_;
  }
  return *this;
}

static inline int MatrixMemError() {
  mprinterr("Error: Not enough memory to allocate pair-wise matrix.\n"
            "Error: Consider using the 'sieve' keyword to reduce memory usage.\n");
  return 1;
}

/** Intended for use with cluster pairwise distance calculations
  * where frames may be sieved. The underlying TriangleMatrix will
  * only be set up to hold the actual number of frames based on
  * the sieve value. The Ignore array will be set up based on
  * the original number of frames.
  */
int DataSet_Cmatrix::SetupWithSieve(size_t sizeIn, size_t sieveIn, int iseed)
{
  if (sievedFrames_.SetSieve( sieveIn, sizeIn, iseed )) return 1;
  // Sieved distances should be ignored.
  if (sievedFrames_.Type() != ClusterSieve::NONE) {
    // Set up the ignore array to ignore sieved frames
    ignore_.assign(sizeIn, true);
    size_t actual_nrows = 0;
    for (size_t frame = 0; frame < sizeIn; frame++)
      if (sievedFrames_.FrameToIdx(frame) != -1) {
        ignore_[frame] = false;
        ++actual_nrows;
      }
    // Set up underlying TriangleMatrix for sieved frames.
    mprintf("\tEstimated pair-wise matrix memory usage: > %s\n",
            ByteString(Mat_.sizeInBytes( 0L, actual_nrows ), BYTE_DECIMAL).c_str());
    try { Mat_.resize( 0L, actual_nrows ); }
    catch (const std::bad_alloc&) { return MatrixMemError(); }
    mprintf("\tPair-wise matrix set up with sieve, %zu frames, %zu sieved frames.\n",
            sizeIn, actual_nrows);
  } else {
    mprintf("\tEstimated pair-wise matrix memory usage: > %s\n",
            ByteString(Mat_.sizeInBytes( 0L, sizeIn ), BYTE_DECIMAL).c_str());
    try { Mat_.resize( 0L, sizeIn ); }
    catch (const std::bad_alloc&) { return MatrixMemError(); }
    ignore_.assign(sizeIn, false);
    mprintf("\tPair-wise matrix set up, %zu frames\n", sizeIn);
  }
  return 0;
}

// DataSet_Cmatrix::SetupIgnore()
int DataSet_Cmatrix::SetupIgnore(size_t originalNrows, std::vector<char> const& ignoreIn,
                                 int sieveIn)
{
  ignore_.assign( originalNrows, false );
  if (!ignoreIn.empty()) {
    if (originalNrows != ignoreIn.size()) {
      mprinterr("Internal Error: Original # rows %zu != ignore array size %zu\n",
                originalNrows, ignoreIn.size());
      return 1;
    }
    for (size_t row = 0; row < originalNrows; ++row)
      if (ignoreIn[row] == 'T')
        ignore_[row] = true;
  }
  // Setup sieve class
  if (sievedFrames_.SetSieve( sieveIn, ignore_ )) {
    mprinterr("Error: Could not set sieve from DataSet_Cmatrix file.\n");
    return 1;
  }
  mprintf("\tSet up %s: %u original rows, %u actual rows, %u elements, sieve=%i\n",
          legend(), originalNrows, Mat_.Nrows(), Mat_.size(), sieveIn);
  return 0;
}

// DataSet_Cmatrix::SetupMatrix()
int DataSet_Cmatrix::SetupMatrix(size_t sizeIn) {
  if (Mat_.resize( 0L, sizeIn )) return 1;
  ignore_.assign( sizeIn, false );
  //sieve_ = 1;
  return 0;
}

// DataSet_Cmatrix::FindMin()
/** Find the minimum; set corresponding row and column. Cannot currently
  * be used for sieved frames.
  */
double DataSet_Cmatrix::FindMin(int& iOut, int& jOut) const {
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

// DataSet_Cmatrix::PrintElements()
void DataSet_Cmatrix::PrintElements() const {
  if (sievedFrames_.MaxFrames()==0) {
    // This is ClusterDistances matrix. Ignore size is # of sieved frames, sievedFrames is 0.
    unsigned int iVal = 0;
    unsigned int jVal = 1;
    for (size_t idx = 0UL; idx < Nelements(); ++idx) {
      if (!ignore_[iVal] && !ignore_[jVal])
        mprintf("\t%u %u %f\n",iVal,jVal,Mat_[idx]);
      // Increment indices
      jVal++;
      if (jVal >= ignore_.size()) {
        iVal++;
        jVal = iVal + 1;
      }
    }
  } else {
    // This is FrameDistances matrix. Ignore and sievedFrames size is # of original frames.
    for (unsigned int row = 0; row != Nframes(); row++)
      for (unsigned int col = row + 1; col != Nframes(); col++)
        if (!ignore_[row] && !ignore_[col])
          mprintf("\t%u %u %f\n", row+1, col+1, GetFdist(col, row));
  }
}

// DataSet_Cmatrix::DataSize()
size_t DataSet_Cmatrix::DataSize() const {
  return ( Mat_.DataSize() +
           (ignore_.capacity()*sizeof(bool) + sizeof(ignore_)) +
           sievedFrames_.DataSize() );
}

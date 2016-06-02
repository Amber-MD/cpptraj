#include "DataSet_Cmatrix_MEM.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // ByteString

void DataSet_Cmatrix_MEM::WriteBuffer(CpptrajFile& outfile, SizeArray const& pIn) const {
  size_t x = (size_t)pIn[0];
  size_t y = (size_t)pIn[1];
  if ( x >= Mat_.Ncols() || y >= Mat_.Nrows() )
    outfile.Printf(format_.fmt(), 0.0);
  else
    outfile.Printf(format_.fmt(), Mat_.element(x,y));
}

int DataSet_Cmatrix_MEM::Allocate(SizeArray const& sizeIn) {
  int err = 0;
  if (!sizeIn.empty()) {
    // Sanity check.
    if (sizeIn.size() > 1 && sizeIn[1] != sizeIn[0])
      mprintf("Warning: DataSet_Cmatrix dimensions must be equal (%zu != %zu)\n"
              "Warning: Matrix will be %zu x %zu upper triangle\n",
              sizeIn[0], sizeIn[1], sizeIn[0], sizeIn[0]);
    err = Mat_.resize(0, sizeIn[0]); // Upper triangle
  } else
    Mat_.clear();
  return err;
}

// -----------------------------------------------------------------------------
/*void DataSet_Cmatrix_MEM::Clear() {
  Mat_.clear();
  ignore_.clear();
  sievedFrames_.Clear();
}*/

// DataSet_Cmatrix_MEM::AllocateCmatrix()
int DataSet_Cmatrix_MEM::AllocateCmatrix(size_t sizeIn)
{
  // Set up underlying TriangleMatrix for given number of frames.
  mprintf("\tEstimated pair-wise matrix memory usage: > %s\n",
          ByteString(Mat_.sizeInBytes( 0L, sizeIn ), BYTE_DECIMAL).c_str());
  try { Mat_.resize( 0L, sizeIn ); }
  catch (const std::bad_alloc&) {
    mprinterr("Error: Not enough memory to allocate pair-wise matrix.\n"
              "Error: Consider using the 'sieve' keyword to reduce memory usage.\n");
    return 1;
  }
  return 0;
}

// DataSet_Cmatrix::SetupIgnore()
/*
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
*/

// DataSet_Cmatrix_MEM::DataSize()
size_t DataSet_Cmatrix_MEM::DataSize() const {
  return ( Mat_.DataSize() +
//           (ignore_.capacity()*sizeof(bool) + sizeof(ignore_)) +
           sievedFrames_.DataSize() );
}

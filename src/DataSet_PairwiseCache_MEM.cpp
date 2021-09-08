#include "DataSet_PairwiseCache_MEM.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

int DataSet_PairwiseCache_MEM::Allocate(SizeArray const& sizeIn) {
  int err = 0;
  if (!sizeIn.empty()) {
    // Sanity check.
    if (sizeIn.size() > 1 && sizeIn[1] != sizeIn[0])
      mprintf("Warning: DataSet_PairwiseCache dimensions must be equal (%zu != %zu)\n"
              "Warning: Matrix will be %zu x %zu upper triangle\n",
              sizeIn[0], sizeIn[1], sizeIn[0], sizeIn[0]);
    err = Mat_.resize(0, sizeIn[0]); // Upper triangle
  } else
    Mat_.clear();
  return err;
}


int DataSet_PairwiseCache_MEM::SetupCache(unsigned int Ntotal, Cframes const& framesToCache,
                                          int sieve, std::string const& metricDescription)
{
  unsigned int sizeIn = framesToCache.size();
  if (sizeIn > 0) {
    mprintf("\tEstimated pair-wise matrix memory usage: >= %s\n",
            ByteString(Mat_.sizeInBytes( 0L, sizeIn ), BYTE_DECIMAL).c_str());
    try { Mat_.resize(0, sizeIn); } // Upper triangle
    catch (const std::bad_alloc&) {
      mprinterr("Error: Not enough memory to allocate pair-wise matrix.\n"
                "Error: Consider using the 'sieve' keyword to reduce memory usage.\n");
      return 1;
    }
#   ifdef DEBUG_CLUSTER
    mprintf("DEBUG: PairwiseMatrix_MEM set up for %i rows, size= %zu bytes.\n",
            Mat_.Nrows(), Mat_.sizeInBytes());
#   endif
  } else
    Mat_.clear();
  SetSieveVal( sieve );
  SetMetricDescrip( metricDescription );
  return SetupFrameToIdx(framesToCache, Ntotal);
}

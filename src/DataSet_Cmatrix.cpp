#include "DataSet_Cmatrix.h"
#include "CpptrajStdio.h"

void DataSet_Cmatrix::PrintElements() const {
  // NOTE: Matrix is square upper triangle, Nrows == Ncols
  for (unsigned int row = 0; row != Nrows(); row++)
    for (unsigned int col = row + 1; col != Nrows(); col++)
      mprintf("\t%u %u %f\n", row+1, col+1, GetFdist(col, row));
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
  // Allocate based on actual size of the matrix.
  if (AllocateCmatrix( sievedFrames_.ActualNframes() )) return 1;
  if (sievedFrames_.Type() != ClusterSieve::NONE)
    mprintf("\tPair-wise matrix set up with sieve, %zu frames, %i sieved frames.\n",
            sievedFrames_.MaxFrames(), sievedFrames_.ActualNframes());
  else
    mprintf("\tPair-wise matrix set up, %zu frames\n", sizeIn);
  return 0;
}

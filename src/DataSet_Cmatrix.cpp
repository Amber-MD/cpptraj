#include "DataSet_Cmatrix.h"
#include "CpptrajStdio.h"

void DataSet_Cmatrix::PrintElements() const {
  // NOTE: Matrix is square upper triangle, Nrows == Ncols
  for (unsigned int row = 0; row != Nrows(); row++)
    for (unsigned int col = row + 1; col != Nrows(); col++)
      mprintf("\t%u %u %f\n", row+1, col+1, GetFdist(col, row));
}

/** Set up sieving info as necessary and set up cluster based on actual
  * number of frames to be clustered.
  */
int DataSet_Cmatrix::SetupWithSieve(ClusterDist* CdistIn, size_t sizeIn, int sieveIn, int iseed)
{
  if (CdistIn == 0) {
    mprinterr("Internal Error: DataSet_Cmatrix::SetupWithSieve called with empty ClusterDist.\n");
    return 1;
  }
  metricDescription_.assign( CdistIn->Description() );
  if (sievedFrames_.SetSieve( sieveIn, sizeIn, iseed )) return 1;
  if (AllocateCmatrix( sievedFrames_.ActualNframes() )) return 1;
  if (SetCdist(CdistIn)) return 1;
  if (sievedFrames_.Type() != ClusterSieve::NONE)
    mprintf("\tPair-wise matrix set up with sieve, %zu frames, %i sieved frames.\n",
            sievedFrames_.MaxFrames(), sievedFrames_.ActualNframes());
  else
    mprintf("\tPair-wise matrix set up, %zu frames\n", sizeIn);
  return 0;
}

/** Set up sieve info from an array that contains 'T' if the frame was sieved
  * out and 'F' otherwise.
  */
int DataSet_Cmatrix::SetSieveFromArray(std::vector<char> const& sieveStatus, int sieveIn)
{
  if (sieveStatus.empty()) return 1;
  // Setup sieve class
  if (sievedFrames_.SetSieve( sieveIn, sieveStatus )) {
    mprinterr("Error: Could not set sieve from cluster matrix file.\n");
    return 1;
  }
  mprintf("\tSet up %s: %zu original frames, %u actual frames, %zu elements",
          legend(), sievedFrames_.MaxFrames(), sievedFrames_.ActualNframes(), Nelements());
  if (sievedFrames_.Type() == ClusterSieve::REGULAR)
    mprintf(", sieve= %i.\n", sievedFrames_.Sieve());
  else if (sievedFrames_.Type() == ClusterSieve::RANDOM)
    mprintf(", random sieve.\n");
  else
    mprintf(".\n");
  return 0;
}

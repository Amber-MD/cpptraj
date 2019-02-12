#include "DataSet_PairwiseCache_NC.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // ByteString

int DataSet_PairwiseCache_NC::SetupCache(unsigned int Ntotal, Cframes const& framesToCache,
                                         int sieve, std::string const& metricDescription)
{
  if (Meta().Fname().empty()) {
    mprinterr("Internal Error: NetCDF pairwise cache file name not set.\n");
    return 1;
  }
  unsigned int sizeIn = framesToCache.size();
  mprintf("\tPairwise cache file: '%s'\n", Meta().Fname().full());
  mprintf("\tEstimated pair-wise matrix disk usage: > %s\n",
          ByteString( ((sizeIn*(sizeIn-1))/2)*sizeof(float), BYTE_DECIMAL).c_str());
    if (file_.CreateCmatrix(Meta().Fname(), Ntotal, sizeIn, sieve, metricDescription))
      return 1;
  // Write actual frames array if necessary
  if (sieve != 1) {
    if (file_.WriteFramesArray( framesToCache ))
      return 1;
  }
  // Reopen in SHARE mode for random access
  if (file_.ReopenSharedWrite( Meta().Fname() )) return 1;
  // TODO - these are saved in the NetCDF file as well - remove redundancy?
  SetSieveVal( sieve );
  SetMetricDescrip( metricDescription );

  return SetupFrameToIdx(framesToCache, Ntotal);
}

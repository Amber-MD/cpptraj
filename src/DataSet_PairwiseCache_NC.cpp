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

  return SetupFrameToIdx(framesToCache, Ntotal);
}


/** Print cached distances to stdout. */
void DataSet_PairwiseCache_NC::PrintCached() const {
  for (Cframes::const_iterator it1 = FrameToIdx().begin(); it1 != FrameToIdx().end(); ++it1)
  {
    if (*it1 != -1) {
      for (Cframes::const_iterator it2 = it1 + 1; it2 != FrameToIdx().end(); ++it2)
      {
        if (*it2 != -1)
          mprintf("\t%i %i %f\n", *it1+1, *it2+1, file_.GetCmatrixElement(*it1, *it2));
      }
    }
  }
}

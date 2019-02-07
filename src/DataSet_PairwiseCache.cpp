#include "DataSet_PairwiseCache.h"
#include "CpptrajStdio.h"

using namespace Cpptraj;
using namespace Cluster;

DataSet_PairwiseCache::DataSet_PairwiseCache() :
  sieve_(0)
{}

DataSet_PairwiseCache::DataSet_PairwiseCache(DataType t) :
  DataSet(t, PWCACHE, TextFormat(TextFormat::DOUBLE, 12, 4), 2),
  sieve_(0)
{}


/** Set up frame number to matrix index for caching. This should be called
  * by each inheriting DataSet_PairwiseCache's SetupCache routine.
  */
int DataSet_PairwiseCache::SetupFrameToIdx(Cframes const& framesToCache, unsigned int Ntotal)
{
  frameToIdx_.assign(Ntotal, -1);
  int idx = 0;
  for (Cframes_it it = framesToCache.begin(); it != framesToCache.end(); ++it)
    frameToIdx_[*it] = idx++;
# ifdef DEBUG_CLUSTER
  // DEBUG
  mprintf("DEBUG: frameToMat\n");
  for (Cframes_it it = frameToIdx_.begin(); it != frameToIdx_.end(); ++it)
    mprintf("\tframeToIdx_[%u] = %i\n", it - frameToIdx_.begin(), *it);
# endif
  return 0;
}

/** Check that the given frames match the already-cached frames. */
bool DataSet_PairwiseCache::CachedFramesMatch(Cframes const& frames) const
{
  Cframes_it it0 = frames.begin();
  for (Cframes_it it1 = frameToIdx_.begin(); it0 != frameToIdx_.end(); ++it1)
  {
    if (*it1 != -1) {
      if (it0 == frames.end() || *it1 != *it0) return false;
      ++it0;
    }
  }
  return true;
}

/** In given StatusArray, T means frame is present, F means not present. */
int DataSet_PairwiseCache::SetupFromStatus(StatusArray const& frameIsPresent, int sieveIn)
{
  frameToIdx_.clear();
  frameToIdx_.reserve( frameIsPresent.size() );
  int idx = 0;
  for (StatusArray::const_iterator it = frameIsPresent.begin();
                                   it != frameIsPresent.end(); ++it)
    if (*it == 'T')
      frameToIdx_.push_back( idx++ );
    else
      frameToIdx_.push_back( -1 );
  sieve_ = sieveIn;
  return 0;
}

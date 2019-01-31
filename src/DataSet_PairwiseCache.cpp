#include "DataSet_PairwiseCache.h"
#include "CpptrajStdio.h"

using namespace Cpptraj;
using namespace Cluster;

/** Set up frame number to matrix index for caching. */
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


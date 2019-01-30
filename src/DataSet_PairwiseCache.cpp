#include "DataSet_PairwiseCache.h"
#include "CpptrajStdio.h"

using namespace Cpptraj;
using namespace Cluster;

/** Set up frame number to matrix index for caching. */
int DataSet_PairwiseCache::setupFrameToMat(Cframes const& framesToCache, unsigned int Ntotal)
{
  frameToMat_.assign(Ntotal, -1);
  int idx = 0;
  for (Cframes_it it = framesToCache.begin(); it != framesToCache.end(); ++it)
    frameToMat_[*it] = idx++;
# ifdef DEBUG_CLUSTER
  // DEBUG
  mprintf("DEBUG: frameToMat\n");
  for (Cframes_it it = frameToMat_.begin(); it != frameToMat_.end(); ++it)
    mprintf("\tframeToMat_[%u] = %i\n", it - frameToMat_.begin(), *it);
# endif
  return 0;
}


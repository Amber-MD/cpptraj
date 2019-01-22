#include "PairwiseMatrix_MEM.h"
#ifdef DEBUG_CLUSTER
#include "../CpptrajStdio.h" // DEBUG
#endif

/** \return distance between frames (cached or uncached). */
double Cpptraj::Cluster::PairwiseMatrix_MEM::Frame_Distance(int f1, int f2) const {
  int idx1 = frameToMat_[f1];
  if (idx1 != -1) {
    int idx2 = frameToMat_[f2];
    if (idx2 != -1) {
      return Mat_.element(idx1, idx2);
    }
  }
  // If here, distance was not cached.
  return MetricPtr()->FrameDist(f1, f2);
}

/** Requests that distances between given frames be cached in memory. */
int Cpptraj::Cluster::PairwiseMatrix_MEM::CacheDistances(Cframes const& framesToCache) {
  Mat_.resize(0, framesToCache.size());
# ifdef DEBUG_CLUSTER
  mprintf("DEBUG: PairwiseMatrix_MEM set up for %i rows, size= %zu bytes.\n",
          Mat_.Nrows(), Mat_.sizeInBytes());
# endif
  frameToMat_.assign(DistMetric().Ntotal(), -1);
  int idx = 0;
  for (Cframes_it it = framesToCache.begin(); it != framesToCache.end(); ++it)
    frameToMat_[*it] = idx++;
# ifdef DEBUG_CLUSTER
  // DEBUG
  mprintf("DEBUG: frameToMat\n");
  for (Cframes_it it = frameToMat_.begin(); it != frameToMat_.end(); ++it)
    mprintf("\tframeToMat_[%u] = %i\n", it - frameToMat_.begin(), *it);
# endif
  return CalcFrameDistances( framesToCache );
}

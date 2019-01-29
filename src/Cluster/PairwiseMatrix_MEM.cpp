#include "PairwiseMatrix_MEM.h"
#include "../CpptrajStdio.h"

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
  if (setupFrameToMat( framesToCache )) return 1;

  return CalcFrameDistances( framesToCache );
}

/** Print cached distances to stdout. */
void Cpptraj::Cluster::PairwiseMatrix_MEM::PrintCached() const {
  for (Cframes::const_iterator it1 = frameToMat_.begin(); it1 != frameToMat_.end(); ++it1)
  {
    if (*it1 != -1) {
      for (Cframes::const_iterator it2 = it1 + 1; it2 != frameToMat_.end(); ++it2)
      {
        if (*it2 != -1)
          mprintf("\t%i %i %f\n", *it1+1, *it2+1, Mat_.element(*it1, *it2));
      }
    }
  }
}


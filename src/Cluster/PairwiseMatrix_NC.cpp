#include "PairwiseMatrix_NC.h"
#include "../CpptrajStdio.h"
#include "../StringRoutines.h" // ByteString

/** \return distance between frames (cached or uncached). */
double Cpptraj::Cluster::PairwiseMatrix_NC::Frame_Distance(int f1, int f2) const {
  int idx1 = frameToMat_[f1];
  if (idx1 != -1) {
    int idx2 = frameToMat_[f2];
    if (idx2 != -1) {
      return  file_.GetCmatrixElement( idx1, idx2 );
    }
  }
  // If here, distance was not cached.
  return MetricPtr()->FrameDist(f1, f2);
}

/** Requests that distances between given frames be cached in memory. */
int Cpptraj::Cluster::PairwiseMatrix_NC::CacheDistances(Cframes const& framesToCache) {
  if (fname_.empty()) {
    mprinterr("Internal Error: PairwiseMatrix_NC::CacheDistances(): File name not set.\n");
    return 1;
  }
  size_t sizeIn = framesToCache.size();
  mprintf("\tPairwise cache file: '%s'\n", fname_.full());
  mprintf("\tEstimated pair-wise matrix disk usage: > %s\n",
          ByteString( ((sizeIn*(sizeIn-1))/2)*sizeof(float), BYTE_DECIMAL).c_str());
  if (file_.CreateCmatrix(fname_, DistMetric().Ntotal(), sizeIn, // TODO sieve value
                          1, DistMetric().Description()))
    return 1;
  // Write actual frames array if necessary // TODO enable
  //if (sievedFrames_.Type() != ClusterSieve::NONE) {
  //  if (file_.WriteFramesArray( framesToCache ))
  //    return 1;
  //}
  // Reopen in SHARE mode for random access
  if (file_.ReopenSharedWrite( fname_ )) return 1;

  if (setupFrameToMat( framesToCache )) return 1;

  return CalcFrameDistances( framesToCache );
}

#include "PairwiseMatrix.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Metric.h"
#include "../CpptrajStdio.h"
#include "../DataSet_PairwiseCache.h"
#include "../ProgressBar.h"

/** Set up PairwiseMatrix with Metric and optional cache. */
int Cpptraj::Cluster::PairwiseMatrix::Setup(Metric* metric, DataSet_PairwiseCache* cache)
{
  if (metric == 0) return 1;
  metric_ = metric;
  cache_  = cache;
  return 0;
}

/** \return distance between frames (cached or uncached). */ // TODO inline?
double Cpptraj::Cluster::PairwiseMatrix::Frame_Distance(int f1, int f2) const {
  if (cache_ != 0)
  {
    // TODO protect against f1/f2 out of bounds.
    int idx1 = cache_->FrameToIdx()[f1];
    if (idx1 != -1) {
      int idx2 = cache_->FrameToIdx()[f2];
      if (idx2 != -1)
        return cache_->CachedDistance(idx1, idx2);
    }
  }
  // If here, distance was not cached or no cache.
  return metric_->FrameDist(f1, f2);
}

/** Request that distances for frames in given array be cached.
  * \param framesToCache the frames to cache.
  * \param sieveIn Sieve value (if any) used to generate framesToCache. This is
  *        purely for bookkeeping inside DataSet_PairwiseCache.
  * \param mismatch_fatal When true, if a cache is present but the frames to
  *        cache do not match what is in the cache, exit with an error;
  *        otherwise recalculate the cache.
  */
int Cpptraj::Cluster::PairwiseMatrix::CacheDistances(Cframes const& framesToCache, int sieveIn,
                                                     bool mismatch_fatal)
{
  if (framesToCache.size() < 1) return 0;
  // If no cache we can leave.
  if (cache_ == 0) return 0;

  bool do_cache = true;
  if (cache_->Size() > 0) {
    do_cache = false;
    mprintf("\tUsing existing cache '%s'\n", cache_->legend());
    // If cache is already populated, check that it is valid.
    // The frames to cache must match cached frames.
    if (!cache_->CachedFramesMatch( framesToCache )) {
      if (mismatch_fatal) {
        mprinterr("Error: Frames to cache do not match those in existing cache.\n");
        return 1;
      } else {
        mprintf("Warning: Frames to cache do not match those in existing cache.\n"
                "Warning: Re-calculating pairwise cache.\n");
        do_cache = true;
      }
    }
    // TODO Check metric? Total frames?
  }

  if (do_cache) {
    // Sanity check
    if (metric_ == 0) {
      mprinterr("Internal Error: PairwiseMatrix::CacheDistances(): Metric is null.\n");
      return 1;
    }
    // Cache setup
    if (cache_->SetupCache( metric_->Ntotal(), framesToCache, sieveIn, metric_->Description() ))
      return 1;
    // Fill cache
    if (CalcFrameDistances(framesToCache))
      return 1;
  }
  return 0;
}

/** Cache distances between given frames using SetElement(). */
int Cpptraj::Cluster::PairwiseMatrix::CalcFrameDistances(Cframes const& framesToCache)
{
  mprintf("\tCaching distances for %zu frames.\n", framesToCache.size());

  int f2end = (int)framesToCache.size();
  int f1end = f2end - 1;
  ParallelProgress progress(f1end);
  int f1, f2;
  // For OMP, every other thread will need its own Cdist.
  Metric* MyMetric = metric_;
# ifdef _OPENMP
# pragma omp parallel private(MyMetric, f1, f2) firstprivate(progress)
  {
  int mythread = omp_get_thread_num();
  progress.SetThread( mythread );
  if (mythread == 0) {
    mprintf("\tParallelizing pairwise distance calc with %i threads\n", omp_get_num_threads());
    MyMetric = metric_;
  } else
    MyMetric = metric_->Copy();
# pragma omp for schedule(dynamic)
# endif
  for (f1 = 0; f1 < f1end; f1++) {
    progress.Update(f1);
    for (f2 = f1 + 1; f2 < f2end; f2++)
      cache_->SetElement( f1, f2, MyMetric->FrameDist(framesToCache[f1], framesToCache[f2]) );
  }
# ifdef _OPENMP
  if (mythread > 0)
    delete MyCdist;
  } // END omp parallel
# endif
  progress.Finish();
  return 0;
}

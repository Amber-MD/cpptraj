#include "PairwiseMatrix.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "../ProgressBar.h"
#include "../CpptrajStdio.h"

/*
Cpptraj::Cluster::PairwiseMatrix::PairwiseMatrix(Type t, Metric* metric) :
  type_(t),
  metric_(metric)
{}*/

/** Set up frame number to matrix index for caching. */
/*
int Cpptraj::Cluster::PairwiseMatrix::setupFrameToMat(Cframes const& framesToCache)
{
  if (metric_ == 0) {
    mprinterr("Internal Error: PairwiseMatrix::setupFrameToMat(): null metric.\n");
    return 1;
  }
  frameToMat_.assign(metric_->Ntotal(), -1);
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
*/

int Cpptraj::Cluster::PairwiseMatrix::CacheDistances(Cframes const& framesToCache) {
  if (framesToCache.size() < 1) return 0;
  if (cache_ == 0) {
    mprinterr("Internal Error: PairwiseMatrix::CacheDistances(): Cache set is null.\n");
    return 1;
  }
  if (metric_ == 0) {
    mprinterr("Internal Error: PairwiseMatrix::CacheDistances(): Metric is null.\n");
    return 1;
  }
  return CalcFrameDistances(framesToCache);
}

/** Cache distances between given frames using SetElement(). */
int Cpptraj::Cluster::PairwiseMatrix::CalcFrameDistances(Cframes const& framesToCache)
{
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


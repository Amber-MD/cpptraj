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

const char DataSet_PairwiseCache::PRESENT_ = 'F';

const char DataSet_PairwiseCache::ABSENT_ = 'T';

/** Set up frame number to matrix index for caching. This should be called
  * by each inheriting DataSet_PairwiseCache's SetupCache routine.
  */
int DataSet_PairwiseCache::SetupFrameToIdx(Cframes const& framesToCache, unsigned int Ntotal)
{
  frameToIdx_.assign(Ntotal, -1);
  int idx = 0;
  for (Cframes_it it = framesToCache.begin(); it != framesToCache.end(); ++it)
    frameToIdx_[*it] = idx++;
//# ifdef DEBUG_CLUSTER
  // DEBUG
  mprintf("DEBUG: frameToMat\n");
  for (Cframes_it it = frameToIdx_.begin(); it != frameToIdx_.end(); ++it)
    mprintf("\tframeToIdx_[%u] = %i\n", it - frameToIdx_.begin(), *it);
//# endif
  return 0;
}

/** Check that the given frames match the already-cached frames. */
bool DataSet_PairwiseCache::CachedFramesMatch(Cframes const& framesIn) const
{
  Cframes_it frm0 = framesIn.begin();
  for (int frm1 = 0; frm1 != (int)frameToIdx_.size(); frm1++)
  {
    if (frameToIdx_[frm1] != -1) {
      if (frm0 == framesIn.end() || *frm0 != frm1) return false;
      ++frm0;
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
    if (*it == PRESENT_)
      frameToIdx_.push_back( idx++ );
    else
      frameToIdx_.push_back( -1 );
  sieve_ = sieveIn;
  return 0;
}

/** \return Array containing frames present in the pairwise cache. */
DataSet_PairwiseCache::Cframes DataSet_PairwiseCache::PresentFrames() const {
  Cframes presentFrames;
  presentFrames.reserve( Nrows() );
  int frm = 0;
  for (frm = 0; frm != (int)frameToIdx_.size(); ++frm)
    if (frameToIdx_[frm] != -1)
      presentFrames.push_back( frm );
  return presentFrames;
}

/** Print cached distances to stdout. */
void DataSet_PairwiseCache::PrintCached() const {
  for (Cframes::const_iterator it1 = FrameToIdx().begin(); it1 != FrameToIdx().end(); ++it1)
  {
    if (*it1 != -1) {
      for (Cframes::const_iterator it2 = it1 + 1; it2 != FrameToIdx().end(); ++it2)
      {
        if (*it2 != -1)
          mprintf("\t%zu %i %f\n", it1-FrameToIdx().begin()+1, it2-FrameToIdx().begin()+1,
                  CachedDistance(*it1, *it2));
      }
    }
  }
}

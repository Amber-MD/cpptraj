#include "Sieve.h"
#include "../DataSet_PairwiseCache.h"
#include "../Random.h"

/** CONSTRUCTOR */
Cpptraj::Cluster::Sieve::Sieve() :
  type_(NONE),
  sieve_(1)
{}

/** Determine the sieving type based on the given sieve value. */
void Cpptraj::Cluster::Sieve::DetermineTypeFromSieve( int sieveIn ) {
  sieve_ = sieveIn;
  // Determine sieve type from sieve value.
  if (sieve_ < -1)
    type_ = RANDOM;
  else if (sieve_ < 2) {
    type_ = NONE;
    sieve_ = 1;
  } else
    type_ = REGULAR;
}

/** Setup which frames should be clustered.*/
int Cpptraj::Cluster::Sieve::SetFramesToCluster(int sieveIn, std::size_t maxFrames, int iseed)
{
  // Sanity check. Should never be called with maxFrames < 1
  if (maxFrames < 1) return 1;
  DetermineTypeFromSieve( sieveIn );
  frameIsPresent_.clear();
  framesToCluster_.clear();
  sievedOut_.clear();
  // ---------------------------------------------
  if (type_ == NONE) 
  {
    // No sieving; frame == index
    framesToCluster_.reserve( maxFrames );
    for (unsigned int i = 0; i < maxFrames; i++)
      framesToCluster_.push_back( i );
  }
  // ---------------------------------------------
  else if (type_ == REGULAR)
  {
    // Regular sieveing; index = frame / sieve
    framesToCluster_.reserve( maxFrames / sieve_ + 1 );
    sievedOut_.reserve( maxFrames - (maxFrames / sieve_) );
    unsigned int tgt = 0;
    for (unsigned int i = 0; i < maxFrames; i++)
    {
      if (i == tgt) {
        framesToCluster_.push_back( i );
        tgt += sieve_;
      } else
        sievedOut_.push_back( i );
    }
  }
  // ---------------------------------------------
  else if (type_ == RANDOM)
  {
    // Random sieving; maxframes / sieve random indices
    Cframes frameToIdx( maxFrames, -1 );
    framesToCluster_.reserve( maxFrames / sieve_ + 1 );
    sievedOut_.reserve( maxFrames - (maxFrames / sieve_) );
    unsigned int maxidx = maxFrames - 1;
    Random_Number random;
    random.rn_set( iseed );
    for (unsigned int i = 0; i < maxFrames; i -= sieve_)
    {
      bool frame_generated = false;
      // Pick until we pick a frame that has not yet been selected.
      while (!frame_generated) {
        unsigned int iframe = random.rn_num_interval(0, maxidx);
        if (frameToIdx[iframe] == -1) {
          frameToIdx[iframe] = 1;
          frame_generated = true;
        }
      }
    }
    // Put indices in order
    for (unsigned int i = 0; i < maxFrames; i++)
      if (frameToIdx[i] == 1)
        framesToCluster_.push_back( i );
      else
        sievedOut_.push_back( i );
  }
  // ---------------------------------------------
//  MakeIdxToFrame();
  return 0;
}

/** Set frames to cluster and sieved out frames from pairwise cache. */
int Cpptraj::Cluster::Sieve::SetupFromCache(DataSet_PairwiseCache const& cache,
                                            std::size_t maxFrames)
{
  if (cache.Size() < 1) {
    //mprinterr("Error: Cannot setup frames to cluster from empty cache.\n");
    return 1;
  }

  frameIsPresent_.clear();
  framesToCluster_.clear();
  sievedOut_.clear();
  DetermineTypeFromSieve( cache.SieveVal() );

  if (type_ == NONE && (cache.Nrows() != maxFrames)) {
    // # of frames in cache does not match frames to cluster
    return 2;
  }

  unsigned int frm = 0;
  for (; frm != cache.FrameToIdx().size(); frm++)
  {
    if (cache.FrameToIdx()[frm] == -1)
      sievedOut_.push_back( frm );
    else
      framesToCluster_.push_back( frm );
  }
  // Check that # of sieved frames in the cache matches # to cluster
  if (cache.Nrows() != framesToCluster_.size()) {
    //mprinterr("Error: # frames in cache (%zu) != # frames to cluster (%zu).\n",
    //          cache.Nrows(), maxFrames);
    framesToCluster_.clear();
    sievedOut_.clear();
    return 2;
  }
  // Anything left is consiedered sieved out.
  for (; frm < maxFrames; frm++)
    sievedOut_.push_back( frm );
  return 0;
}

/** Clear the sieve. */
void Cpptraj::Cluster::Sieve::Clear() {
  framesToCluster_.clear();
  sievedOut_.clear();
  frameIsPresent_.clear();
  type_ = NONE;
  sieve_ = 1;
}

/** Create an array containing 'true' for frames that are present, 'false' otherwise.
  */
void Cpptraj::Cluster::Sieve::GenerateFrameIsPresentArray() {
  if (frameIsPresent_.empty()) {
    // Start everything out as false TODO search for max val?
    frameIsPresent_.assign( framesToCluster_.size() + sievedOut_.size(), false );
    // Set true for frames to cluster
    for (Cframes::const_iterator it = framesToCluster_.begin();
                                 it != framesToCluster_.end(); ++it)
      frameIsPresent_[*it] = true;
  }
}

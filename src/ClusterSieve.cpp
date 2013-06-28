#include "ClusterSieve.h"
#include "Random.h"

// CONSTRUCTOR
ClusterSieve::ClusterSieve() : type_(NONE), sieve_(1) {}

// ClusterSieve::SetSieve()
int ClusterSieve::SetSieve(int sieveIn, size_t maxFrames, int iseed) {
  if (maxFrames < 1) return 1;
  sieve_ = sieveIn;
  // Determine sieve type from sieve value.
  if (sieve_ == -1)
    type_ = RANDOM;
  else if (sieve_ < 2) {
    type_ = NONE;
    sieve_ = 1;
  } else
    type_ = REGULAR;
  frameToIdx_.clear();
  if (type_ == NONE) 
  { // No sieving; frame == index
    frameToIdx_.reserve( maxFrames );
    for (unsigned int i = 0; i < maxFrames; i++)
      frameToIdx_.push_back( i );
  }
  else if (type_ == REGULAR)
  { // Regular sieveing; index = frame / sieve
    frameToIdx_.assign( maxFrames, -1 );
    int idx = 0;
    for (unsigned int i = 0; i < maxFrames; i += sieve_)
      frameToIdx_[i] = idx++;
  }
  else if (type_ == RANDOM)
  { // Random sieving; maxframes / sieve random indices
    frameToIdx_.assign( maxFrames, -1 );
    double dmax = (double)maxFrames;
    Random_Number random;
    random.rn_set( iseed );
    int idx = 0;
    for (unsigned int i = 0; i < maxFrames; i += sieve_)
    {
      bool frame_generated = false;
      // Pick until we pick a frame that has not yet been selected.
      while (!frame_generated) {
        double dframe = dmax * random.rn_gen();
        int iframe = (int)dframe;
        if (frameToIdx_[iframe] == -1) {
          frameToIdx_[iframe] = idx++;
          frame_generated = true;
        }
      }
    }
    // TODO: Put frames in order?
  }
  return 0;
}

// ClusterSieve::Frames()
ClusterSieve::SievedFrames ClusterSieve::Frames() const {
  SievedFrames frames;
  for (unsigned int frame = 0; frame != frameToIdx_.size(); frame++)
  {
    if (frameToIdx_[frame] != -1)
      frames.push_back( frame );
  }
  return frames;
}  

#include <algorithm>
#include "Cframes.h"
#include "../Random.h"

bool Cpptraj::Cluster::Cframes::HasFrame(int frame) const {
  Iarray::const_iterator it = std::find(frames_.begin(), frames_.end(), frame);
  return !(it == frames_.end());
}

void Cpptraj::Cluster::Cframes::Remove(int frame) {
  Iarray::iterator pend = std::remove( frames_.begin(), frames_.end(), frame);
  size_t newsize = pend - frames_.begin();
  frames_.resize( newsize );
}

void Cpptraj::Cluster::Cframes::Sort() {
  std::sort(frames_.begin(), frames_.end());
}

void Cpptraj::Cluster::Cframes::DetermineTypeFromSieve( int sieveIn ) {
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
int Cpptraj::Cluster::Cframes::SetFramesToCluster(int sieveIn, size_t maxFrames, int iseed)
{
  if (maxFrames < 1) return 1;
  DetermineTypeFromSieve( sieveIn );
//  frameToIdx_.clear();
  // ---------------------------------------------
  if (type_ == NONE) 
  { // No sieving; frame == index
    frames_.reserve( maxFrames );
    for (unsigned int i = 0; i < maxFrames; i++)
      frames_.push_back( i );
    //actualNframes_ = (int)maxFrames;
  }
  // ---------------------------------------------
  else if (type_ == REGULAR)
  { // Regular sieveing; index = frame / sieve
    //frames_.assign( maxFrames, -1 );
    frames_.reserve( maxFrames / sieve_ + 1 );
    //int idx = 0;
    for (unsigned int i = 0; i < maxFrames; i += sieve_)
      frames_.push_back( i );
      //frameToIdx_[i] = idx++;
    //actualNframes_ = idx;
  }
  // ---------------------------------------------
  else if (type_ == RANDOM)
  { // Random sieving; maxframes / sieve random indices
    Iarray frameToIdx( maxFrames, -1 );
    frames_.reserve( maxFrames / sieve_ + 1 );
    double dmax = (double)maxFrames;
    Random_Number random;
    random.rn_set( iseed );
    for (unsigned int i = 0; i < maxFrames; i -= sieve_)
    {
      bool frame_generated = false;
      // Pick until we pick a frame that has not yet been selected.
      while (!frame_generated) {
        double dframe = dmax * random.rn_gen();
        int iframe = (int)dframe;
        if (frameToIdx[iframe] == -1) {
          frameToIdx[iframe] = 1;
          frame_generated = true;
        }
      }
    }
    // Put indices in order
    //int idx = 0;
    for (unsigned int i = 0; i < maxFrames; i++)
      if (frameToIdx[i] == 1)
        frames_.push_back( i );
    //actualNframes_ = idx;
  }
  // ---------------------------------------------
//  MakeIdxToFrame();
  return 0;
}


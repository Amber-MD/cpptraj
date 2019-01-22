#include "Sieve.h"
#include "../Random.h"

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
int Cpptraj::Cluster::Sieve::SetFramesToCluster(int sieveIn, size_t maxFrames, int iseed)
{
  if (maxFrames < 1) return 1;
  DetermineTypeFromSieve( sieveIn );
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

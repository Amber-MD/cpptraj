#include "ClusterNode.h"

// CONSTRUCTOR
ClusterNode::ClusterNode() :
  avgClusterDist_(0),
  eccentricity_(0),
  num_(0),
  centroid_(0)
{}

// Set initial centroid to front, even though that will probably be wrong
// when number of frames in the list > 1
ClusterNode::ClusterNode(std::list<int>& frameListIn, int numIn) :
  avgClusterDist_(0),
  eccentricity_(0),
  num_(numIn),
  centroid_(frameListIn.front()),
  frameList_(frameListIn)
{}

ClusterNode::ClusterNode(const ClusterNode& rhs) :
  avgClusterDist_( rhs.avgClusterDist_ ),
  eccentricity_( rhs.eccentricity_ ),
  num_( rhs.num_ ),
  centroid_( rhs.centroid_ ),
  frameList_( rhs.frameList_ )
{}

ClusterNode& ClusterNode::operator=(const ClusterNode& rhs) {
  if (&rhs == this) return *this;
  avgClusterDist_ = rhs.avgClusterDist_;
  eccentricity_ = rhs.eccentricity_;
  num_ = rhs.num_;
  centroid_ = rhs.centroid_;
  frameList_ = rhs.frameList_;
  return *this;
}

/// Used to sort clusters by # of frames in cluster
/** Use > since we give higher priority to larger clusters. */
bool ClusterNode::operator<(const ClusterNode& rhs) const {
  return ( frameList_.size() > rhs.frameList_.size() );
}

void ClusterNode::MergeFrames( ClusterNode& rhs) {
  frameList_.splice( frameList_.begin(), rhs.frameList_ );
}

void ClusterNode::SetNum(int numIn) {
  num_ = numIn;
  frameList_.sort();
}

/** When sieving, frame #s are off by sieve. Fix frame and centroid #s. */
void ClusterNode::FrameSieveOffset(int sieve) {
  for (std::list<int>::iterator fnum = frameList_.begin();
                                fnum != frameList_.end(); ++fnum)
    (*fnum) *= sieve;
  centroid_ *= sieve;
}

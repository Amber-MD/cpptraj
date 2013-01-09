#include <cfloat> // DBL_MAX
#include "ClusterNode.h"

// CONSTRUCTOR
ClusterNode::ClusterNode() :
  avgClusterDist_(0),
  eccentricity_(0),
  num_(0),
  centroidframe_(0),
  centroid_(0)
{}

// DESTRUCTOR
ClusterNode::~ClusterNode() {
  if (centroid_ != 0) delete centroid_;
}

// Set initial centroid frame to front, even though that will probably be wrong
// when number of frames in the list > 1
ClusterNode::ClusterNode(ClusterDist* Cdist, std::list<int> const& frameListIn, int numIn) :
  avgClusterDist_(0.0),
  eccentricity_(0.0),
  num_(numIn),
  centroidframe_(frameListIn.front()),
  frameList_(frameListIn),
  centroid_(Cdist->NewCentroid(frameList_))
{}

// COPY CONSTRUCTOR
ClusterNode::ClusterNode(const ClusterNode& rhs) :
  avgClusterDist_( rhs.avgClusterDist_ ),
  eccentricity_( rhs.eccentricity_ ),
  num_( rhs.num_ ),
  centroidframe_( rhs.centroidframe_ ),
  frameList_( rhs.frameList_ ),
  centroid_(0)
{
  if (rhs.centroid_ != 0)
    centroid_ = rhs.centroid_->Copy();
}

// ASSIGNMENT
ClusterNode& ClusterNode::operator=(const ClusterNode& rhs) {
  if (&rhs == this) return *this;
  avgClusterDist_ = rhs.avgClusterDist_;
  eccentricity_ = rhs.eccentricity_;
  num_ = rhs.num_;
  centroidframe_ = rhs.centroidframe_;
  frameList_ = rhs.frameList_;
  if (centroid_ != 0) delete centroid_;
  if (rhs.centroid_ != 0)
    centroid_ = rhs.centroid_->Copy();
  else
    centroid_ = 0;
  return *this;
}

/** Use > since we give higher priority to larger clusters. */
bool ClusterNode::operator<(const ClusterNode& rhs) const {
  return ( frameList_.size() > rhs.frameList_.size() );
}

/** Frames from rhs go to this cluster. rhs frames are removed. */
void ClusterNode::MergeFrames( ClusterNode& rhs) {
  frameList_.splice( frameList_.begin(), rhs.frameList_ );
}

/** Set internal cluster number, sort frame list. */
void ClusterNode::SetNum(int numIn) {
  num_ = numIn;
  frameList_.sort();
}

/** Find the frame in the given cluster that is the centroid, i.e. has the
  * lowest cumulative distance to every other point in the cluster.
  */
int ClusterNode::FindCentroidFrame(ClusterMatrix const& FrameDistancesIn) {
  double mindist = DBL_MAX;
  int minframe = -1;
  for (frame_iterator frm1 = frameList_.begin(); frm1 != frameList_.end(); ++frm1)
  {
    double cdist = 0.0;
    for (frame_iterator frm2 = frameList_.begin(); frm2 != frameList_.end(); ++frm2)
    {
      if (frm1 == frm2) continue;
      cdist += FrameDistancesIn.GetElement(*frm1, *frm2);
    }
    if (cdist < mindist) {
      mindist = cdist;
      minframe = *frm1;
    }
  }
  if (minframe == -1) 
    return 1;
  centroidframe_ = minframe;
  return 0;
}

/** Calculate the eccentricity of this cluster (i.e. the largest distance
  * between any two points in the cluster).
  */
void ClusterNode::CalcEccentricity(ClusterMatrix const& FrameDistancesIn) {
  double maxdist = 0.0;
  frame_iterator frame1_end = frameList_.end();
  --frame1_end;
  for (frame_iterator frm1 = frameList_.begin(); frm1 != frameList_.end(); ++frm1)
  {
    frame_iterator frm2 = frm1;
    ++frm2;
    for (; frm2 != frameList_.end(); ++frm2) {
      double fdist = FrameDistancesIn.GetElement(*frm1, *frm2);
      if (fdist > maxdist)
        maxdist = fdist;
    }
  }
  eccentricity_ = maxdist;
}

/** Calculate average distance between all members in cluster and
  * the centroid. 
  */
double ClusterNode::CalcAvgToCentroid( ClusterDist* Cdist )
{
  double avgdist = 0.0;
  //int idx = 0; // DEBUG
  //mprintf("AVG DISTANCES FOR CLUSTER %d:\n", Num()); // DEBUG
  for (frame_iterator frm = frameList_.begin(); frm != frameList_.end(); ++frm)
  {
    double dist = Cdist->FrameCentroidDist( *frm, centroid_ );
    //mprintf("\tDist to %i is %f\n", idx++, dist); // DEBUG
    avgdist += dist;
  }
  return ( avgdist / (double)frameList_.size() );
}

#include <cfloat> // DBL_MAX
#include "ClusterNode.h"
#include "DataSet_double.h"
#include "DS_Math.h"
// DEBUG
#include "CpptrajStdio.h"

// CONSTRUCTOR
ClusterNode::ClusterNode() :
  avgClusterDist_(0),
  eccentricity_(0),
  num_(0),
  centroid_(0)
{}

// Set initial centroid to front, even though that will probably be wrong
// when number of frames in the list > 1
ClusterNode::ClusterNode(std::list<int> const& frameListIn, int numIn) :
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

/** Compute the centroid (avg) coords for each atom from all frames in this
  * cluster.
  */
// NOTE: Use mask?
void ClusterNode::CalcCentroidFrame(DataSet_Coords const& coordsIn, AtomMask const& maskIn) 
{
  cframe_.SetupFrame( maskIn.Nselected() );
  cframe_.ZeroCoords();
  // Temp frame to hold input coords
  Frame frameIn( maskIn.Nselected() );
  for (frame_iterator frm = frameList_.begin(); frm != frameList_.end(); ++frm)
  {
    frameIn.SetFromCRD( coordsIn[ *frm ], coordsIn.NumBoxCrd(), maskIn );
    cframe_ += frameIn;
  }
  cframe_.Divide( (double)frameList_.size() );
  //mprintf("\t\tFirst 3 centroid coords (of %i): %f %f %f\n", cframe_.Natom(), cframe_[0],
  //        cframe_[1],cframe_[2]);
}

/** Find the frame in the given cluster that is the centroid, i.e. has the
  * lowest cumulative distance to every other point in the cluster.
  */
int ClusterNode::FindCentroid(ClusterMatrix const& FrameDistancesIn) {
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
  centroid_ = minframe;
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

/** Get the average distance between all frames in the cluster */
void ClusterNode::CalcAvgFrameDist(ClusterMatrix const& FrameDistancesIn) {
  DataSet_double distances;
  int numframes = (int)frameList_.size(); 
  distances.Resize( ((numframes * numframes) - numframes) / 2 );
  int nd = 0;
  frame_iterator frame1_end = frameList_.end();
  --frame1_end;
  for (frame_iterator frm1 = frameList_.begin(); frm1 != frameList_.end(); ++frm1)
  {
    frame_iterator frm2 = frm1;
    ++frm2;
    for (; frm2 != frameList_.end(); ++frm2)
      distances[nd++] = FrameDistancesIn.GetElement( *frm1, *frm2 );
  }
  internalAvg_ = DS_Math::Avg( distances, &internalSD_ );
}

/** Calculate average distance between all frames in cluster and
  * the centroid frame.
  */
// TODO: usemass, DME, etc
double ClusterNode::CalcAvgToCentroid( DataSet_Coords const& coordsIn, AtomMask const& maskIn ) 
{
  // TODO: Check that mask size matches centroid
  // Temp frame to hold input coords
  Frame frameIn( maskIn.Nselected() );
  // Temp frame to hold reference coords
  Frame c_ref( cframe_.Natom() );
  double avgdist = 0.0;
  //int idx = 0; // DEBUG
  for (frame_iterator frm = frameList_.begin(); frm != frameList_.end(); ++frm)
  {
    frameIn.SetFromCRD( coordsIn[ *frm ], coordsIn.NumBoxCrd(), maskIn );
    c_ref.SetCoordinates( cframe_ );
    double dist = frameIn.RMSD( c_ref, false ); // Best-fit RMSD
    //mprintf("\tDist to %i is %f\n", idx++, dist);
    avgdist += dist;
  }
  return ( avgdist / (double)frameList_.size() );
}

double ClusterNode::CentroidDist( ClusterNode const& rhs ) {
  // Temp frame for this centroid
  Frame cent1 = cframe_;
  // Temp frame for other centroid
  Frame cent2 = rhs.cframe_;
  // distance
  double dist = cent1.RMSD( cent2, false );
  return dist;
}

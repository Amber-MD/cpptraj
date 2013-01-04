#include <cfloat> // DBL_MAX
#include <cmath> // sqrt
#include "ClusterNode.h"
#include "DataSet_Coords.h"
#include "DataSet_double.h"
#include "DS_Math.h"
// DEBUG
#include "CpptrajStdio.h"

// CONSTRUCTOR
ClusterNode::ClusterNode() :
  avgClusterDist_(0),
  internalAvg_(0.0),
  internalSD_(0.0),
  eccentricity_(0),
  cval_(0.0),
  num_(0),
  centroid_(0)
{}

// Set initial centroid to front, even though that will probably be wrong
// when number of frames in the list > 1
ClusterNode::ClusterNode(std::list<int> const& frameListIn, int numIn) :
  avgClusterDist_(0.0),
  internalAvg_(0.0),
  internalSD_(0.0),
  eccentricity_(0.0),
  cval_(0.0),
  num_(numIn),
  centroid_(frameListIn.front()),
  frameList_(frameListIn)
{}

/// Create singleton cluster whose only purpose will be centroid calc (e.g. in sieving)
ClusterNode::ClusterNode(DataSet* dsIn, int frameIn, RMSoptions const& rmsopt) :
  avgClusterDist_(0.0),
  internalAvg_(0.0),
  internalSD_(0.0),
  eccentricity_(0.0),
  cval_(0.0),
  num_(frameIn),
  centroid_(0)
{
  if (dsIn->Type() == DataSet::COORDS) {
    DataSet_Coords* coords = (DataSet_Coords*)dsIn;
    cframe_.SetupFrame( rmsopt.mask.Nselected() );
    coords->GetFrame( frameIn, cframe_, rmsopt.mask );
  } else
    cval_ = dsIn->Dval( frameIn );
}

// COPY CONSTRUCTOR
ClusterNode::ClusterNode(const ClusterNode& rhs) :
  avgClusterDist_( rhs.avgClusterDist_ ),
  internalAvg_( rhs.internalAvg_ ),
  internalSD_( rhs.internalSD_ ),
  eccentricity_( rhs.eccentricity_ ),
  cval_( rhs.cval_ ),
  num_( rhs.num_ ),
  centroid_( rhs.centroid_ ),
  frameList_( rhs.frameList_ ),
  cframe_( rhs.cframe_ )
{}

// ASSIGNMENT
ClusterNode& ClusterNode::operator=(const ClusterNode& rhs) {
  if (&rhs == this) return *this;
  avgClusterDist_ = rhs.avgClusterDist_;
  internalAvg_ = rhs.internalAvg_;
  internalSD_ = rhs.internalSD_;
  eccentricity_ = rhs.eccentricity_;
  num_ = rhs.num_;
  centroid_ = rhs.centroid_;
  frameList_ = rhs.frameList_;
  cframe_ = rhs.cframe_;
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

/** When sieving, frame #s are off by sieve. Fix frame and centroid #s. */
void ClusterNode::FrameSieveOffset(int sieve) {
  for (std::list<int>::iterator fnum = frameList_.begin();
                                fnum != frameList_.end(); ++fnum)
    (*fnum) *= sieve;
  centroid_ *= sieve;
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

/** Two potential modes based on type of DataSet. If DataSet is COORDS, 
  * Compute the centroid (avg) coords for each atom from all frames in this
  * cluster. In order to have more representative avg coords, RMS fit to  
  * centroid as it is being built.
  * If other type of dataset just calc average of double vals.
  */
// NOTE: Should ONLY be fitting if RMS was previously calcd with FIT!
void ClusterNode::CalculateCentroid(DataSet* dsIn, RMSoptions const& rmsopt) 
{
  if (dsIn->Type() == DataSet::COORDS) {
    DataSet_Coords* coords = (DataSet_Coords*)dsIn;
    Matrix_3x3 Rot;
    Vec3 Trans;
    // Reset atom count for centroid.
    cframe_.ClearAtoms();
    // Frame to hold input coords
    Frame frameIn( rmsopt.mask.Nselected() );
    for (frame_iterator frm = frameList_.begin(); frm != frameList_.end(); ++frm)
    {
      coords->GetFrame( *frm, frameIn, rmsopt.mask );
      if (cframe_.empty()) {
        cframe_ = frameIn;
        if (!rmsopt.nofit)
          cframe_.CenterOnOrigin(rmsopt.useMass);
      } else {
        if (!rmsopt.nofit) {
          frameIn.RMSD_CenteredRef( cframe_, Rot, Trans, rmsopt.useMass );
          frameIn.Rotate( Rot );
        }
        cframe_ += frameIn;
      }
    }
    cframe_.Divide( (double)frameList_.size() );
    //mprintf("\t\tFirst 3 centroid coords (of %i): %f %f %f\n", cframe_.Natom(), cframe_[0],
    //        cframe_[1],cframe_[2]);
  } else {
    // Generic DataSet; get average of all frames; empty centroid frame
    cframe_.ClearAtoms();
    cval_ = 0.0;
    for (frame_iterator frm = frameList_.begin(); frm != frameList_.end(); ++frm)
      cval_ += dsIn->Dval( *frm );
    cval_ /= (double)frameList_.size();
  }
}


/** Calculate average distance between all members in cluster and
  * the centroid. 
  */
double ClusterNode::CalcAvgToCentroid( DataSet* dsIn, RMSoptions const& rmsopt) 
{
  if (dsIn->Type() == DataSet::COORDS) {
    DataSet_Coords* coords = (DataSet_Coords*)dsIn;
    double dist;
    // TODO: Check that mask size matches centroid
    // Temp frame to hold input coords
    Frame frameIn( rmsopt.mask.Nselected() );
    double avgdist = 0.0;
    int idx = 0; // DEBUG
    for (frame_iterator frm = frameList_.begin(); frm != frameList_.end(); ++frm)
    {
      coords->GetFrame( *frm, frameIn, rmsopt.mask);
      if (rmsopt.useDME)
        dist = frameIn.DISTRMSD( cframe_ );
      else {
        if (!rmsopt.nofit)
          // Centroid is already at origin.
          dist = frameIn.RMSD_CenteredRef( cframe_, rmsopt.useMass ); // Best-fit RMSD
        else
          dist = frameIn.RMSD_NoFit( cframe_, rmsopt.useMass );
      }
      mprintf("\tDist to %i is %f\n", idx++, dist); // DEBUG
      avgdist += dist;
    }
    return ( avgdist / (double)frameList_.size() );
  } else {
    // Generic DataSet
    double sumdiff2 = 0.0;
    for (frame_iterator frm = frameList_.begin(); frm != frameList_.end(); ++frm)
    {
      double diff = dsIn->Dval( *frm ) - cval_;
      sumdiff2 += (diff * diff);
    }
    sumdiff2 /= (double)frameList_.size();
    return sqrt( sumdiff2 );
  }
}

/** Calculate distance between centroid frames. */
double ClusterNode::CentroidDist( ClusterNode const& rhs, RMSoptions const& rmsopt ) {
  if (!cframe_.empty()) {
    if (rmsopt.useDME)
      return cframe_.DISTRMSD( rhs.cframe_ );
    else {
      if (!rmsopt.nofit)
        // Centroids are already at origin.
        return cframe_.RMSD_CenteredRef( rhs.cframe_, rmsopt.useMass );
      else
        return cframe_.RMSD_NoFit( rhs.cframe_, rmsopt.useMass );
    }
  } else {
    double diff = cval_ - rhs.cval_;
    return fabs( diff );
  }
}

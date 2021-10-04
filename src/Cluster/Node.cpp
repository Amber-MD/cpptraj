//#include <cfloat> // DBL_MAX
#include "Node.h"
#include "MetricArray.h"
#include "../DataSet_float.h"
#include "../DataSet_integer.h"

// CONSTRUCTOR
Cpptraj::Cluster::Node::Node() :
  avgSil_(0),
  eccentricity_(0),
  refRms_(0),
  num_(-1),
  needsUpdate_(true)
{}

// DESTRUCTOR
Cpptraj::Cluster::Node::~Node() { }

/** Create new cluster with given number containing given frames. Calculate
  * initial centroid and set initial best rep frame to front, even though 
  * that will probably be wrong when number of frames in the list > 1.
  */
Cpptraj::Cluster::Node::Node(MetricArray& Cdist, Cframes const& frameListIn, int numIn) :
  frameList_(frameListIn),
  bestReps_(1, RepPair(frameListIn.front(), 0.0)),
  eccentricity_(0.0),
  num_(numIn),
  needsUpdate_(true)
{
  Cdist.NewCentroid( centroids_, frameListIn );
}

// COPY CONSTRUCTOR
Cpptraj::Cluster::Node::Node(const Node& rhs) :
  frameList_( rhs.frameList_ ),
  centroids_(rhs.centroids_),
  bestReps_( rhs.bestReps_ ),
  eccentricity_( rhs.eccentricity_ ),
  num_( rhs.num_ ),
  needsUpdate_( rhs.needsUpdate_ )
{ }

// ASSIGNMENT
Cpptraj::Cluster::Node& Cpptraj::Cluster::Node::operator=(const Node& rhs) {
  if (&rhs == this) return *this;
  eccentricity_ = rhs.eccentricity_;
  num_ = rhs.num_;
  bestReps_ = rhs.bestReps_;
  frameList_ = rhs.frameList_;
  centroids_ = rhs.centroids_;
  needsUpdate_ = rhs.needsUpdate_;
  return *this;
}

/** Calculate centroid of frames in this Node using given metric. */
void Cpptraj::Cluster::Node::CalculateCentroid(MetricArray& Cdist) {
  if (centroids_.empty())
    Cdist.NewCentroid( centroids_, frameList_ );
  else
    Cdist.CalculateCentroid( centroids_, frameList_ );
}

/** Calculate the eccentricity of this cluster (i.e. the largest distance
  * between any two points in the cluster).
  */
void Cpptraj::Cluster::Node::CalcEccentricity(MetricArray& FrameDistancesIn) {
  double maxdist = 0.0;
  frame_iterator frame1_end = frameList_.end();
  --frame1_end;
  for (frame_iterator frm1 = frameList_.begin(); frm1 != frameList_.end(); ++frm1)
  {
    frame_iterator frm2 = frm1;
    ++frm2;
    for (; frm2 != frameList_.end(); ++frm2) {
      double fdist = FrameDistancesIn.Frame_Distance(*frm1, *frm2);
      if (fdist > maxdist)
        maxdist = fdist;
    }
  }
  eccentricity_ = maxdist;
}

/** Calculate average distance between all members in cluster and
  * the centroid. 
  */
double Cpptraj::Cluster::Node::CalcAvgToCentroid( MetricArray& Cdist ) const
{
  double avgdist = 0.0;
  //int idx = 0; // DEBUG
  //mprintf("AVG DISTANCES FOR CLUSTER %d:\n", Num()); // DEBUG
  for (frame_iterator frm = frameList_.begin(); frm != frameList_.end(); ++frm)
  {
    double dist = Cdist.FrameCentroidDist( *frm, centroids_ );
    //mprintf("\tDist to %i is %f\n", idx++, dist); // DEBUG
    avgdist += dist;
  }
  return ( avgdist / (double)frameList_.size() );
}

/** Remove specified frame, update the centroid. */
void Cpptraj::Cluster::Node::RemoveFrameUpdateCentroid(MetricArray& Cdist, int frame) {
  Cdist.FrameOpCentroid(frame, centroids_, (double)frameList_.size(),
                        Metric::SUBTRACTFRAME);
  RemoveFrameFromCluster( frame );
}

/** Add specified frame, update the centroid. */
void Cpptraj::Cluster::Node::AddFrameUpdateCentroid(MetricArray& Cdist, int frame) {
  Cdist.FrameOpCentroid(frame, centroids_, (double)frameList_.size(),
                        Metric::ADDFRAME);
  AddFrameToCluster( frame );
}

/** Calculate cluster population vs time. */
void Cpptraj::Cluster::Node::CalcCpopVsTime(DataSet_float& pvt, unsigned int maxFrames,
                                            CnormType normType)
const
{
  // TODO clear pvt array?
  float pop = 0.0;
  // Loop over all frames in cluster
  for (frame_iterator f = beginframe(); f != endframe(); ++f)
  {
    if (*f > (int)pvt.Size())
      pvt.Resize( *f, pop );
    pop = pop + 1.0;
    pvt[*f] = pop;
  }
  // Ensure pop v time set is maxFrames long
  if (pvt.Size() < maxFrames)
    pvt.Resize( maxFrames, pop );
  // Normalization
  if (normType == CLUSTERPOP) {
    float norm = 1.0 / (float)Nframes();
    for (unsigned int frm = 0; frm < maxFrames; ++frm)
      pvt[frm] = pvt[frm] * norm;
  } else if (normType == FRAME) {
    float norm = 1.0;
    for (unsigned int frm = 0; frm < maxFrames; ++frm)
    {
      pvt[frm] = pvt[frm] / norm;
      norm = norm + 1.0;
    }
  }
}

/** Create cluster "lifetime" set, with 1 for cluster present and 0 for absent. */
void Cpptraj::Cluster::Node::CreateLifetimeSet(DataSet_integer& life, unsigned int maxFrames)
const
{
  life.Resize( maxFrames );
  for (frame_iterator f = beginframe(); f != endframe(); ++f)
    life.SetElement( *f, 1);
}

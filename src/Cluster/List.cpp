#include <limits>
#include "List.h"
#include "../CpptrajStdio.h"
#include "../ProgressBar.h"

void Cpptraj::Cluster::List::PrintClusters() const {
  //mprintf("CLUSTER: %u clusters, %u frames.\n", clusters_.size(),
  //        FrameDistances().OriginalNframes() );
  mprintf("CLUSTER: %u clusters.\n", clusters_.size());
  for (cluster_iterator C = begincluster(); C != endcluster(); C++) {
    mprintf("\t%8i : ",C->Num());
    for (Node::frame_iterator fnum = C->beginframe();
                              fnum != C->endframe(); ++fnum)
      mprintf("%i,",(*fnum)+1);
    mprintf("\n");
  }
}

/** Create cluster number vs time data set. */
int Cpptraj::Cluster::List::CreateCnumVsTime(DataSet_integer* ds, unsigned int maxFrames)
const
{
  if (ds == 0) {
    mprinterr("Internal Error: CreateCnumVsTime() called with null data set\n");
    return 1;
  }
  DataSet_integer& cnum_temp = static_cast<DataSet_integer&>( *ds );
  cnum_temp.Resize( maxFrames );
  // Make all clusters start at -1. This way cluster algorithms that
  // have noise points (i.e. no cluster assigned) will be distinguished.
  std::fill(cnum_temp.begin(), cnum_temp.end(), -1);

  for (cluster_iterator C = begincluster(); C != endcluster(); C++)
  {
    //mprinterr("Cluster %i:\n",CList->CurrentNum());
    int cnum = C->Num();
    // Loop over all frames in the cluster
    for (Node::frame_iterator frame = C->beginframe(); frame != C->endframe(); frame++)
    {
      //mprinterr("%i,",*frame);
      cnum_temp[ *frame ] = cnum;
    }
    //mprinterr("\n");
    //break;
  }
  return 0;
}

/** Sort clusters by population and renumber. */
int Cpptraj::Cluster::List::Sort() {
  clusters_.sort();
  int newNum = 0;
  for (cluster_it it = clusters_.begin(); it != clusters_.end(); ++it)
    it->SetNum( newNum++ );
  return 0;
}

/** Update centroids. */
void Cpptraj::Cluster::List::UpdateCentroids( Metric* metric ) {
  for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node) {
    node->SortFrameList();
    node->CalculateCentroid( metric );
  }
}

// -----------------------------------------------------------------------------
void Cpptraj::Cluster::List::AddFramesByCentroid(Cframes const& framesIn, Metric* metricIn)
{
  // NOTE: All cluster centroids must be up to date.
  int idx;
  int nframes = (int)framesIn.size();
  double mindist, dist;
  cluster_it minNode, Cnode;
  ParallelProgress progress( nframes );
  // For OMP, every other thread will need its own Cdist.
  Metric* MyCdist = metricIn;
# ifdef _OPENMP
  // For OMP need a temp. array to hold which frame goes to which cluster to avoid clashes
  std::vector<cluster_it> idxToCluster( nframes, clusters_.end() );
# pragma omp parallel private(MyCdist, idx, dist, mindist, minNode, Cnode) firstprivate(progress)
  {
  int mythread = omp_get_thread_num();
  progress.SetThread( mythread );
  if (mythread == 0) {
    mprintf("\tParallelizing sieve restore calc with %i threads\n", omp_get_num_threads());
    MyCdist = metricIn;
  } else
    MyCdist = metricIn->Copy();
# pragma omp for schedule(dynamic)
# endif
  for (idx = 0; idx < nframes; ++idx) {
    progress.Update( idx );
    int frame = framesIn[idx];
    // Which clusters centroid is closest to this frame?
    mindist = std::numeric_limits<double>::max();
    minNode = clusters_.end();
    for (Cnode = clusters_.begin(); Cnode != clusters_.end(); ++Cnode) {
      dist = MyCdist->FrameCentroidDist(frame, Cnode->Cent());
      if (dist < mindist) {
        mindist = dist;
        minNode = Cnode;
      }
    }
    // Add sieved frame to the closest cluster.
#   ifdef _OPENMP
    idxToCluster[idx] = minNode;
#   else
    minNode->AddFrameToCluster( frame );
#   endif
  } // END loop over frames
# ifdef _OPENMP
  if (mythread > 0)
    delete MyCdist;
  } // END pragma omp parallel
  // Now actually add sieved frames to their appropriate clusters
  for (idx = 0; idx < nframes; idx++)
    idxToCluster[idx]->AddFrameToCluster( framesIn[idx] );
# endif
  progress.Finish();
}


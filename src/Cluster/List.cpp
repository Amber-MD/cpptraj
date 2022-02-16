#include <limits> // double max
#include <list>
#include "List.h"
#include "MetricArray.h"
#include "Node.h"
#include "../CpptrajStdio.h"
#include "../DataSet_integer.h"
#include "../ProgressBar.h"
#ifdef _OPENMP
# include <omp.h>
#endif

void Cpptraj::Cluster::List::PrintClusters() const {
  //mprintf("CLUSTER: %u clusters, %u frames.\n", clusters_.size(),
  //        FrameDistances().OriginalNframes() );
  mprintf("CLUSTER: %zu clusters.\n", clusters_.size());
  for (cluster_iterator C = begincluster(); C != endcluster(); C++) {
    mprintf("\t%8i : ",C->Num());
    for (Node::frame_iterator fnum = C->beginframe();
                              fnum != C->endframe(); ++fnum)
      mprintf("%i,",(*fnum)+1);
    mprintf("\n");
  }
}

/** Create cluster number vs time data set. 
  * \param ds The output integer DataSet
  * \param maxFrames The maximum number of frames.
  * \param offset If specified, add offset to cluster numbers.
  * \param maxCluster If specified, number clusters beyond maxCluster equal to maxCluster
  */
int Cpptraj::Cluster::List::CreateCnumVsTime(DataSet_integer& cnum_temp, unsigned int maxFrames,
                                             int offset, int maxCluster)
const
{
  // Make all clusters start at -1. This way cluster algorithms that
  // have noise points (i.e. no cluster assigned) will be distinguished.
  cnum_temp.Assign( maxFrames, NOISE + offset );

  if (offset == 0 && maxCluster < 0) {
    for (cluster_iterator C = begincluster(); C != endcluster(); C++)
    {
      //mprinterr("Cluster %i:\n",CList->CurrentNum());
      int cnum = C->Num();
      // Loop over all frames in the cluster
      for (Node::frame_iterator frame = C->beginframe(); frame != C->endframe(); frame++)
      {
        //mprinterr("%i,",*frame);
        cnum_temp.SetElement( *frame, cnum );
      }
      //mprinterr("\n");
      //break;
    }
  } else {
    for (cluster_iterator C = begincluster(); C != endcluster(); C++)
    {
      int cnum = C->Num() + offset;
      if (cnum > maxCluster)
        cnum = maxCluster;
       // Loop over all frames in the cluster
      for (Node::frame_iterator frame = C->beginframe(); frame != C->endframe(); frame++)
        cnum_temp.SetElement( *frame , cnum );
    }
  }
  return 0;
}

/** Determine how many different clusters are observed within a given time
  * window.
  */
int Cpptraj::Cluster::List::NclustersObserved(DataSet_integer& clustersVtime,
                                              unsigned int maxFrames,
                                              int windowSize)
const
{
  if (windowSize < 1) {
    mprinterr("Error: Invalid window size for number clusters observed vs time.\n");
    return 1;
  }
  // Determine number of windows
  int nwindows = (int)maxFrames / windowSize;
  if (((int)maxFrames % windowSize) != 0) nwindows++;
  if (debug_ > 0)
    mprintf("DEBUG: %u frames, %i windows, window size %i.\n", maxFrames, nwindows, windowSize);

  // Create a bool array for each window that will record if cluster is present
  // during that window.
  typedef std::vector<bool> Barray;
  typedef std::vector<Barray> Warray;
  Warray Windows;
  Windows.reserve( nwindows );
  for (int w = 0; w != nwindows; w++)
    Windows.push_back( Barray(Nclusters(), false) );

  // Loop over clusters
  for (cluster_iterator C = begincluster(); C != endcluster(); C++)
  {
    // Loop over cluster frames
    for (Node::frame_iterator frame = C->beginframe(); frame != C->endframe(); frame++)
    {
      int wndw = *frame / windowSize;
      Windows[wndw][C->Num()] = true;
    }
  }

  // Count number of unique clusters in each window
  clustersVtime.Allocate( DataSet::SizeArray(1, nwindows) );
  int wndw = 0;
  for (Warray::const_iterator window = Windows.begin(); window != Windows.end(); ++window, ++wndw)
  {
    int nunique = 0;
    for (Barray::const_iterator cpresent = window->begin(); cpresent != window->end(); ++cpresent)
      if (*cpresent) ++nunique;
    clustersVtime.Add(wndw, &nunique);
  }

  clustersVtime.SetDim(Dimension::X, Dimension(windowSize, windowSize, "Frame"));
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
void Cpptraj::Cluster::List::UpdateCentroids( MetricArray& metric ) {
  for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node) {
    node->SortFrameList();
    node->CalculateCentroid( metric );
  }
}

/** Remove clusters with no population. */
void Cpptraj::Cluster::List::RemoveEmptyClusters() {
  cluster_it cnode = clusters_.begin();
  while (cnode != clusters_.end()) {
    if (cnode->Nframes() == 0)
      cnode = clusters_.erase( cnode );
    else
      ++cnode;
  }
}

/** Clear everything. */
void Cpptraj::Cluster::List::Clear() {
  clusters_.clear();
  noise_.clear();
}

// -----------------------------------------------------------------------------
void Cpptraj::Cluster::List::AddFramesByCentroid(Cframes const& framesIn, MetricArray& metricIn)
{
  // NOTE: All cluster centroids must be up to date.
  int idx;
  int nframes = (int)framesIn.size();
  double mindist, dist;
  cluster_it minNode, Cnode;
  ParallelProgress progress( nframes );
  // For OMP, every other thread will need its own Cdist.
  MetricArray* MyCdist = &metricIn;
# ifdef _OPENMP
  // For OMP need a temp. array to hold which frame goes to which cluster to avoid clashes
  std::vector<cluster_it> idxToCluster( nframes, clusters_.end() );
# pragma omp parallel private(MyCdist, idx, dist, mindist, minNode, Cnode) firstprivate(progress)
  {
  int mythread = omp_get_thread_num();
  progress.SetThread( mythread );
  if (mythread == 0) {
    mprintf("\tParallelizing sieve restore calc with %i threads\n", omp_get_num_threads());
    MyCdist = &metricIn;
  } else
    MyCdist = new MetricArray(metricIn);
# pragma omp for
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
    if (idxToCluster[idx] != clusters_.end())
      idxToCluster[idx]->AddFrameToCluster( framesIn[idx] );
# endif
  progress.Finish();
}

/** Add frames to clusters that are within epsilon of a cluster centroid (if
  * sieveToCentroid) or epsilon of any frame in a cluster otherwise.
  */
void Cpptraj::Cluster::List::AddFramesByCentroid(Cframes const& framesIn, MetricArray& metricIn,
                                                 bool sieveToCentroid, double epsilon)
{
  // NOTE: All cluster centroids must be up to date!
  if (sieveToCentroid)
    mprintf("\tRestoring sieved frames if within %.3f of cluster centroid.\n", epsilon);
  else
    mprintf("\tRestoring sieved frames if within %.3f of frame in nearest cluster.\n",
            epsilon);
  // Vars allocated here in case of OpenMP
  int n_sieved_noise = 0;
  int Nsieved = 0;
  int idx;
  int nframes = (int)framesIn.size();
  ParallelProgress progress( nframes );
  // Need a temporary array to hold which frame belongs to which cluster. 
  // Otherwise we could be comparoing sieved frames to other sieved frames.
  std::vector<cluster_it> idxToCluster( nframes, clusters_.end() );
  // For OMP, every other thread will need its own Cdist.
  MetricArray* MyCdist = &metricIn;
# ifdef _OPENMP
# pragma omp parallel private(MyCdist, idx) firstprivate(progress) reduction(+ : Nsieved, n_sieved_noise)
  {
  int mythread = omp_get_thread_num();
  progress.SetThread( mythread );
  if (mythread == 0) {
    mprintf("\tParallelizing calculation with %i threads\n", omp_get_num_threads());
    MyCdist = &metricIn;
  } else
    MyCdist = new MetricArray(metricIn);
# pragma omp for
# endif
  for (idx = 0; idx < nframes; ++idx) {
    progress.Update( idx );
    int frame = framesIn[idx];
    // Which clusters centroid is closest to this frame?
    double mindist = std::numeric_limits<double>::max(); 
    cluster_it minNode = clusters_.end();
    for (cluster_it Cnode = clusters_.begin(); Cnode != clusters_.end(); ++Cnode) {
      double dist = MyCdist->FrameCentroidDist(frame, Cnode->Cent());
      if (dist < mindist) {
        mindist = dist;
        minNode = Cnode;
      }
    }
    bool goodFrame = false;
    if ( mindist < epsilon )
      // Frame is already within epsilon, accept.
      goodFrame = true;
    else if ( !sieveToCentroid ) {
      // Check if any frames in the cluster are closer than epsilon to sieved frame.
      for (int cidx=0; cidx < minNode->Nframes(); cidx++)
      { //TODO just use PairwiseMatrix::Frame_Distance here?
        if ( MyCdist->Frame_Distance(frame, minNode->ClusterFrame(cidx)) < epsilon )
        {
          goodFrame = true;
          break;
        }
      }
    }
    // Add sieved frame to the closest cluster if closest distance is
    // less than epsilon.
    ++Nsieved;
    if ( goodFrame )
      idxToCluster[idx] = minNode;
    else
      n_sieved_noise++;
  } // END loop over frames
# ifdef _OPENMP
  if (mythread > 0)
    delete MyCdist;
  } // END pragma omp parallel
# endif
  progress.Finish();
  // Now actually add sieved frames to their appropriate clusters
  for (idx = 0; idx < nframes; idx++)
    if (idxToCluster[idx] != clusters_.end())
      idxToCluster[idx]->AddFrameToCluster( framesIn[idx] );
  mprintf("\t%i of %i sieved frames were discarded as noise.\n", 
          n_sieved_noise, Nsieved);
}

#include <algorithm> // std::max
#include <limits> // double max
#include "List.h"
#include "MetricArray.h"
#include "Node.h"
#include "../CpptrajStdio.h"
#include "../DataSet_integer.h"
#include "../ProgressBar.h"
#include "../Constants.h" // SMALL TODO use limits?
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

// -----------------------------------------------------------------------------
/** The Davies-Bouldin Index (DBI) measures the average similarity between each
  * cluster and its most similar one; the smaller the DBI, the better. The DBI 
  * is defined as the average, for all clusters X, of fred, where fred(X) = max,
  * across other clusters Y, of (Cx + Cy)/dXY. Here Cx is the average distance
  * from points in X to the centroid, similarly Cy, and dXY is the distance 
  * between cluster centroids.
  * NOTE: To use this, cluster centroids should be fully up-to-date.
  */
double Cpptraj::Cluster::List::ComputeDBI(std::vector<double>& averageDist, MetricArray& metricIn)
const
{
  averageDist.clear();
  averageDist.reserve( clusters_.size() );
  for (cluster_iterator C1 = begincluster(); C1 != endcluster(); ++C1) {
    // Calculate average distance to centroid for this cluster
    averageDist.push_back( C1->CalcAvgToCentroid( metricIn ) );
    //if (outfile.IsOpen())
    //  outfile.Printf("#Cluster %i has average-distance-to-centroid %f\n", 
    //                 C1->Num(), averageDist.back());
  }
  double DBITotal = 0.0;
  unsigned int nc1 = 0;
  for (cluster_iterator c1 = begincluster(); c1 != endcluster(); ++c1, ++nc1) {
    double MaxFred = 0;
    unsigned int nc2 = 0;
    for (cluster_iterator c2 = begincluster(); c2 != endcluster(); ++c2, ++nc2) {
      if (c1 != c2) {
        double Fred = averageDist[nc1] + averageDist[nc2];
        Fred /= metricIn.CentroidDist( c1->Cent(), c2->Cent() );
        if (Fred > MaxFred)
          MaxFred = Fred;
      }
    }
    DBITotal += MaxFred;
  }
  DBITotal /= (double)clusters_.size();
  //if (outfile.IsOpen()) outfile.Printf("#DBI: %f\n", DBITotal);
  return DBITotal;
}

/** The pseudo-F statistic is another measure of clustering goodness. It is 
  * intended to capture the 'tightness' of clusters, and is in essence a ratio
  * of the mean sum of squares between groups to the mean sum of squares within
  * group (Lattin et al., 2003: 291) High values are good. Generally, one 
  * selects a cluster-count that gives a peak in the pseudo-f statistic (or 
  * pSF, for short).
  * Formula: A/B, where A = (T - P)/(G-1), and B = P / (n-G). Here n is the 
  * number of points, G is the number of clusters, T is the total distance from
  * the all-data centroid, and P is the sum (for all clusters) of the distances
  * from the cluster centroid.
  * NOTE: To use this, cluster centroids should be fully up-to-date.
  * NOTE: This calc differs slightly from PTRAJ in that real centroids are used
  *       instead of representative structures.
  */
double Cpptraj::Cluster::List::ComputePseudoF(double& SSRSST, MetricArray& metricIn) const
{
  // Calculation makes no sense with fewer than 2 clusters.
  if (Nclusters() < 2) {
    mprintf("Warning: Fewer than 2 clusters. Not calculating pseudo-F.\n");
    return 0.0;
  }

  // Form a cluster with all points to get a centroid. Use only frames that
  // are in clusters, i.e. ignore noise. Assumes all cluster centroids are
  // up to date.
  Node c_all;
  for (cluster_iterator C1 = begincluster(); C1 != endcluster(); ++C1)
  {
    for (Node::frame_iterator f1 = C1->beginframe(); f1 != C1->endframe(); ++f1)
      c_all.AddFrameToCluster( *f1 );
  }
  // Pseudo-F makes no sense if # clusters == # frames
  if (Nclusters() == c_all.Nframes()) {
    mprintf("Warning: Each frame is in a separate cluster. Not calculating pseudo-F.\n");
    return 0.0;
  }
  c_all.SortFrameList();
  c_all.CalculateCentroid( metricIn );

  // Loop over all clusters
  double gss = 0.0; // between-group sum of squares
  double wss = 0.0; // within-group sum of squares
  for (cluster_iterator C1 = begincluster(); C1 != endcluster(); ++C1)
  {
    for (Node::frame_iterator f1 = C1->beginframe(); f1 != C1->endframe(); ++f1)
    {
      double dist = metricIn.FrameCentroidDist(*f1, c_all.Cent());
      gss += (dist * dist);
      dist = metricIn.FrameCentroidDist(*f1, C1->Cent());
      wss += (dist * dist);
    }
  }
  double d_nclusters = (double)Nclusters();
  double d_ntotal = (double)c_all.Nframes();
  double num = (gss - wss) / (d_nclusters - 1.0);
  double den = wss / (d_ntotal - d_nclusters);
  if (den < Constants::SMALL)
    den = Constants::SMALL;
  double pseudof = num / den;
  if (debug_ > 0)
    mprintf("Pseudo-f: Total distance to centroid is %.4f\n"
            "Pseudo-f: Cluster distance to centroid is %.4f\n"
            "Pseudo-f: Numerator %.4f over denominator %.4f gives %.4f\n", 
            gss, wss, num, den, pseudof);
  //if (outfile.IsOpen()) {
  //  outfile.Printf("#pSF: %f\n", pseudof);
    // This calculation taken directly from ptraj
    SSRSST = pseudof*(d_nclusters-1)/(d_ntotal-d_nclusters+pseudof*(d_nclusters-1));
  //  outfile.Printf("#SSR/SST: %f\n", SSRSST);
  //}

  return pseudof;
}

/** The cluster silhouette is a measure of how well each point fits within
  * a cluster. Values of 1 indicate the point is very similar to other points
  * in the cluster, i.e. it is well-clustered. Values of -1 indicate the point
  * is dissimilar and may fit better in a neighboring cluster. Values of 0
  * indicate the point is on a border between two clusters. 
  */
int Cpptraj::Cluster::List::CalcSilhouette(MetricArray& metrics,
                                            Cframes const& sievedFrames,
                                            bool includeSieved)
{
  for (cluster_it Ci = begin(); Ci != end(); ++Ci)
  {
    Node::SilPairArray& SiVals = Ci->FrameSilhouettes();
    SiVals.clear();
    SiVals.reserve( Ci->Nframes() );
    double avg_si = 0.0;
    int ci_frames = 0;
    for (Node::frame_iterator f1 = Ci->beginframe(); f1 != Ci->endframe(); ++f1)
    {
      if (includeSieved || !sievedFrames.HasFrame( *f1 )) {
        // Calculate the average dissimilarity of this frame with all other
        // points in this frames cluster.
        double ai = 0.0;
        int self_frames = 0;
        if (includeSieved) {
          for (Node::frame_iterator f2 = Ci->beginframe(); f2 != Ci->endframe(); ++f2)
          {
            if (f1 != f2) {
              ai += metrics.Frame_Distance(*f1, *f2);
              ++self_frames;
            }
          }
        } else {
          for (Node::frame_iterator f2 = Ci->beginframe(); f2 != Ci->endframe(); ++f2)
          {
            if (f1 != f2 && !sievedFrames.HasFrame(*f2)) {
              ai += metrics.Frame_Distance(*f1, *f2);
              ++self_frames;
            }
          }
        }
        if (self_frames > 0)
          ai /= (double)self_frames;
        //mprintf("\t\tFrame %i cluster %i ai = %g\n", *f1+1, Ci->Num(), ai);
        // Determine lowest average dissimilarity of this frame with all
        // other clusters.
        double min_bi = std::numeric_limits<double>::max();
        for (cluster_iterator Cj = begincluster(); Cj != endcluster(); ++Cj)
        {
          if (Ci != Cj)
          {
            double bi = 0.0;
            // NOTE: ASSUMING NO EMPTY CLUSTERS
            if (includeSieved) {
              for (Node::frame_iterator f2 = Cj->beginframe(); f2 != Cj->endframe(); ++f2)
                bi += metrics.Frame_Distance(*f1, *f2);
              bi /= (double)Cj->Nframes();
            } else {
              int cj_frames = 0;
              for (Node::frame_iterator f2 = Cj->beginframe(); f2 != Cj->endframe(); ++f2)
              {
                if (!sievedFrames.HasFrame(*f2)) {
                  bi += metrics.Frame_Distance(*f1, *f2);
                  ++cj_frames;
                }
              }
              bi /= (double)cj_frames;
            }
            //mprintf("\t\tFrame %i to cluster %i bi = %g\n", *f1 + 1, Cj->Num(), bi);
            if (bi < min_bi)
              min_bi = bi;
          }
        }
        double max_ai_bi = std::max( ai, min_bi );
        if (max_ai_bi == 0.0)
          mprinterr("Error: Divide by zero in silhouette calculation for frame %i\n", *f1 + 1);
        else {
          double si = (min_bi - ai) / max_ai_bi;
          SiVals.push_back( Node::SilPair(*f1, si) );
          avg_si += si;
          ++ci_frames;
        }
      } // END if frame should be calcd
    } // END loop over cluster frames
    //std::sort( SiVals.begin(), SiVals.end() );
    // DEBUG
    if (debug_ > 1) {
      mprintf("DEBUG: Cluster frame silhouette values for cluster %i\n", Ci->Num());
      for (Node::SilPairArray::const_iterator it = Ci->FrameSilhouettes().begin();
                                              it != Ci->FrameSilhouettes().end(); ++it)
        mprintf("\t%8i %g\n", it->first+1, it->second);
    }
    if (ci_frames > 0)
      avg_si /= (double)ci_frames;
    //mprintf("DEBUG: Cluster silhouette: %8i %g\n", Ci->Num(), avg_si);
    Ci->SetSilhouette( avg_si );
  } // END outer loop over clusters
  return 0;
}

#include <cfloat> // DBL_MAX
#include <vector>
#include "ClusterList.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"
#include "Analysis_Rms2d.h"
#include "ProgressBar.h"

// XMGRACE colors
const char* ClusterList::XMGRACE_COLOR[] = {
  "white", "black", "red", "green", "blue", "yellow", "brown", "grey", "violet",
  "cyan", "magenta", "orange", "indigo", "maroon", "turquoise", "darkgreen"
};

// CONSTRUCTOR
ClusterList::ClusterList() : debug_(0) {}

// ClusterList::SetDebug()
/** Set the debug level */
void ClusterList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0) mprintf("ClusterList debug set to %i\n",debug_);
}

// ClusterList::Renumber()
/** Sort clusters by size and renumber starting from 0, where cluster 0
  * is the largest. Also calculate anything dependent on FrameDistances
  * (i.e. centroid, avg distance within clusters). Renumber frames
  * according to sieve if necessary.
  * NOTE: This destroys indexing into ClusterDistances.
  */
void ClusterList::Renumber(int sieve) {
  // Before clusters are renumbered, calculate the average distance of 
  // this cluster to every other cluster
  double numdist = (double) (clusters_.size() - 1);
  for (cluster_it node = clusters_.begin();
                  node != clusters_.end(); node++)
  {
    double avgclusterdist = 0;
    for (cluster_it node2 = clusters_.begin();
                    node2 != clusters_.end(); node2++)
    {
      if (node == node2) continue;
      //mprintf("DBG:\t\t%i to %i %lf\n",(*node).num, (*node2).num, 
      //        ClusterDistances.GetElement( (*node).num, (*node2).num ));
      avgclusterdist += ClusterDistances_.GetElement( (*node).Num(), 
                                                     (*node2).Num() );
    }
    avgclusterdist /= numdist;
    //mprintf("DBG:\tCluster %i avg dist = %lf\n",(*node).num,avgclusterdist);
    (*node).SetAvgDist( avgclusterdist ); 
  }
  
  clusters_.sort( );
  // Since frame #s may be modifed by sieveing, calculate cluster properties
  // that depend on frame number here.
  int newNum = 0;
  for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node) 
  {
    (*node).SetNum( newNum++ );
    // Find the centroid. Since FindCentroidFrame uses FrameDistances and not
    // ClusterDistances its ok to call after sorting/renumbering.
    if ((*node).FindCentroidFrame( FrameDistances_ )) {
      mprinterr("Error: Could not determine centroid frame for cluster %i\n",
                (*node).Num());
    }
    // Calculate avg and stdev of all distances.
    (*node).CalcAvgFrameDist( FrameDistances_ );
  }

  // If sieveing the frame numbers should actually be multiplied by sieve
  if (sieve > 0) {
    for (cluster_it node = clusters_.begin(); node != clusters_.end(); node++)
      (*node).FrameSieveOffset( sieve );
  }
}

// ClusterList::Summary()
/** Print a summary of clusters.  */
void ClusterList::Summary(std::string const& summaryfile, int maxframesIn) {
  CpptrajFile outfile;
  float fmax = (float)maxframesIn;
  if (outfile.OpenWrite(summaryfile)) {
    mprinterr("Error: ClusterList::Summary: Could not set up file.\n");
    return;
  }

  outfile.Printf("%-8s %8s %8s %8s %8s %8s %8s\n","#Cluster","Frames","Frac",
                     "AvgDist","Stdev","Centroid","AvgCDist");
  for (cluster_it node = clusters_.begin();
                  node != clusters_.end(); node++)
  {
    // Calculate size and fraction of total size of this cluster
    int numframes = (*node).Nframes();
    float frac = (float)numframes / fmax;
    // OUTPUT
    outfile.Printf("%8i %8i %8.3f %8.3lf %8.3lf %8i %8.3lf\n",
                   (*node).Num(), numframes, frac, (*node).InternalAvg(), 
                   (*node).InternalSD(), (*node).Centroid()+1, (*node).AvgDist() );
  } // END loop over clusters
  outfile.CloseFile();
}

// ClusterList::Summary_Half
/** Print a summary of the first half of the data to the second half.
  */
void ClusterList::Summary_Half(std::string const& summaryfile, int maxframesIn) {
  CpptrajFile outfile;
  float fmax = (float)maxframesIn;
  if (outfile.OpenWrite(summaryfile)) {
    mprinterr("Error: ClusterList::Summary_Half: Could not set up file.\n");
    return;
  }

  // Calculate halfway point
  int half = maxframesIn / 2;
  // xmgrace color
  int color = 1;

  outfile.Printf("#%-7s %8s %6s %2s %10s %8s %8s %6s %6s\n", 
                 "Cluster", "Total", "Frac", "C#", "Color", 
                 "NumIn1st", "NumIn2nd","Frac1","Frac2");
  for (cluster_it node = clusters_.begin();
                  node != clusters_.end(); node++)
  {
    // Calculate size and fraction of total size of this cluster
    int numframes = (*node).Nframes();
    float frac = (float)numframes / fmax;
    int numInFirstHalf = 0;
    int numInSecondHalf = 0;
    // DEBUG
    //mprintf("\tCluster %i\n",(*node).num);
    // Count how many frames are in the first half and how many 
    // are in the second half.
    for (ClusterNode::frame_iterator frame1 = (*node).beginframe();
                                     frame1 != (*node).endframe();
                                     frame1++)
    {
      if (*frame1 < half)
        ++numInFirstHalf;
      else
        ++numInSecondHalf;
    }
    float frac1 = (float) numframes;
    frac1 = ((float) numInFirstHalf) / frac1;
    float frac2 = (float) numframes;
    frac2 = ((float) numInSecondHalf) / frac2;
    outfile.Printf("%-8i %8i %6.2f %2i %10s %8i %8i %6.2f %6.2f\n",
                   (*node).Num(), numframes, frac, color, XMGRACE_COLOR[color],
                   numInFirstHalf, numInSecondHalf, frac1, frac2);
    if (color<15) ++color;
  }
  outfile.CloseFile();
}

// -----------------------------------------------------------------------------
// ClusterList::AddCluster()
/** Add a cluster made up of frames specified by the given framelist to 
  * the list.
  */
int ClusterList::AddCluster( std::list<int> const& framelistIn, int numIn ) {
  clusters_.push_back( ClusterNode( framelistIn, numIn ) );
  return 0;
}

// ClusterList::CalcFrameDistances()
int ClusterList::CalcFrameDistances(std::string const& filename, DataSet* dsIn,
                                    DistModeType mode, 
                                    ClusterNode::RMSoptions const& rmsopt)
{
  if (dsIn == 0) {
    mprinterr("Internal Error: ClusterList: Cluster properties DataSet is null.\n");
    return 1;
  }
  if (mode == USE_DATASET)   // Get distances from DataSet
    calcDistFromDataSet( dsIn );
  else if (dsIn->Type() != DataSet::COORDS) {
    mprinterr("Internal Error: DataSet type is not COORDS and mode is not USE_DATASET\n");
    return 1;
  }
  // The following modes only work with COORDS
  if (mode == USE_FILE) {    // Get distances from file
    mprintf(" Loading pair-wise distances from %s\n", filename.c_str());
    if (FrameDistances_.LoadFile( filename, dsIn->Size() )) {
      mprintf("\tLoading pair-wise distances failed - regenerating from frames.\n");
      mode = USE_FRAMES;
    }
  }
  if (mode == USE_FRAMES) { // Get distances from RMSDs between frames.
    if (rmsopt.useDME)
      Analysis_Rms2d::CalcDME( *((DataSet_Coords*)dsIn), FrameDistances_, 
                               rmsopt.mask.MaskString() );
    else
      Analysis_Rms2d::Calc2drms( *((DataSet_Coords*)dsIn), FrameDistances_, 
                                 rmsopt.nofit, rmsopt.useMass, rmsopt.mask.MaskString() );
  }
  
  // Save distances - overwrites old distances
  if (!USE_FILE)
    FrameDistances_.SaveFile( filename );
  // DEBUG - Print Frame distances
  if (debug_ > 1) {
    mprintf("INTIAL FRAME DISTANCES:\n");
    FrameDistances_.PrintElements();
  }
  return 0;
}  
    
// ClusterList::calcDistFromDataSet()
void ClusterList::calcDistFromDataSet( DataSet* ClusterData ) {
  int N = ClusterData->Size();
  //mprintf("DEBUG: xmax is %i\n",N);
  FrameDistances_.Setup(N);
  int max = FrameDistances_.Nelements();
  mprintf(" Calculating distances using dataset %s (%i total).\n",
          ClusterData->Legend().c_str(),max);

  ProgressBar progress(max);
  // LOOP 
  int current = 0;
  for (int i = 0; i < N-1; i++) {
    progress.Update(current);
    double iVal = ClusterData->Dval(i);
    for (int j = i + 1; j < N; j++) {
      double jVal = ClusterData->Dval(j);
      // Calculate abs( delta )
      double delta = iVal - jVal;
      if (delta < 0) delta = -delta;
      FrameDistances_.AddElement( delta );
      current++;
    }
  }
}

// -----------------------------------------------------------------------------
// ClusterList::InitializeClusterDistances()
/** Set up the initial distances between clusters. Should be called before 
  * any clustering is performed. 
  */
void ClusterList::InitializeClusterDistances(LINKAGETYPE linkage) {
  ClusterDistances_.Setup( clusters_.size() );
  // Set up the ignore array
  ClusterDistances_.SetupIgnore();
  // Build initial cluster distances
  if (linkage==AVERAGELINK) {
    for (cluster_it C1_it = clusters_.begin(); 
                    C1_it != clusters_.end(); C1_it++) 
      calcAvgDist(C1_it);
  } else if (linkage==SINGLELINK) {
    for (cluster_it C1_it = clusters_.begin();
                    C1_it != clusters_.end(); C1_it++) 
      calcMinDist(C1_it);
  } else if (linkage==COMPLETELINK) {
    for (cluster_it C1_it = clusters_.begin(); 
                    C1_it != clusters_.end(); C1_it++) 
      calcMaxDist(C1_it);
  }
  if (debug_ > 1) {
    mprintf("CLUSTER: INITIAL CLUSTER DISTANCES:\n");
    ClusterDistances_.PrintElements();
  }
}

// ClusterList::ClusterHierAgglo()
/** Cluster using a hierarchical agglomerative (bottom-up) approach. All frames
  * start in their own cluster. The closest two clusters are merged, and 
  * distances between the newly merged cluster and all remaining clusters are
  * recalculated according to one of the following metrics:
  * - single-linkage  : The minimum distance between frames in clusters are used.
  * - average-linkage : The average distance between frames in clusters are used.
  * - complete-linkage: The maximum distance between frames in clusters are used.
  */
int ClusterList::ClusterHierAgglo(double epsilon, int targetN, LINKAGETYPE linkage) 
{
  mprintf("\tStarting Hierarchical Agglomerative Clustering:\n");
  ProgressBar cluster_progress(-1);
  // Build initial clusters.
  for (int cluster = 0; cluster < FrameDistances_.Nrows(); cluster++) 
    AddCluster( std::list<int>(1, cluster), cluster);
  mprintf("\t%i initial clusters.\n", Nclusters());
  // Build initial cluster distance matrix.
  InitializeClusterDistances(linkage);
  // DEBUG - print initial clusters
  if (debug_ > 1)
    PrintClusters();
  bool clusteringComplete = false;
  int iterations = 0;
  while (!clusteringComplete) {
    // Merge 2 closest clusters. Clustering complete if closest dist > epsilon.
    if (MergeClosest(epsilon, linkage)) break;
    // If the target number of clusters is reached we are done
    if (Nclusters() <= targetN) {
      mprintf("\n\tTarget # of clusters (%i) met (%u), clustering complete.\n", targetN,
              Nclusters());
      break;
    }
    if (Nclusters() == 1) clusteringComplete = true; // Sanity check
    cluster_progress.Update( iterations++ );
  }
  mprintf("\tCompleted after %i iterations, %u clusters.\n",iterations,
          Nclusters());
  return 0;
}

// -----------------------------------------------------------------------------
// ClusterList::PrintClusters()
/** Print list of clusters and frame numbers belonging to each cluster.
  */
void ClusterList::PrintClusters() {
  mprintf("CLUSTER: %u clusters, %i frames.\n", clusters_.size(), FrameDistances_.Nrows() );
  for (cluster_it C = clusters_.begin(); C != clusters_.end(); C++) {
    mprintf("\t%8i : ",(*C).Num());
    for (ClusterNode::frame_iterator fnum = (*C).beginframe();
                                     fnum != (*C).endframe(); ++fnum)
      mprintf("%i,",(*fnum)+1);
    mprintf("\n");
  }
}

// ClusterList::PrintClustersToFile()
/** Print list of clusters in a style similar to ptraj; each cluster is
  * given a line maxframes characters long, with X for each frame that is
  * in the clusters and . for all other frames. Also print out the
  * representative frame numbers.
  */
void ClusterList::PrintClustersToFile(std::string const& filename, int maxframesIn) {
  CpptrajFile outfile;
  std::string buffer;
  
  if ( outfile.OpenWrite(filename) ) {
    mprinterr("Error: PrintClustersToFile: Could not set up file %s\n",
              filename.c_str());
    return;
  }
  outfile.Printf("#Clustering: %u clusters %i frames\n",
                 clusters_.size(), maxframesIn);
  for (cluster_it C1_it = clusters_.begin(); 
                  C1_it != clusters_.end(); C1_it++)
  {
    buffer.clear();
    buffer.resize(maxframesIn, '.');
    for (ClusterNode::frame_iterator frame1 = (*C1_it).beginframe();
                                     frame1 != (*C1_it).endframe();
                                     frame1++)
    {
      buffer[ *frame1 ] = 'X';
    }
    buffer += '\n';
    outfile.Write((void*)buffer.c_str(), buffer.size());
  }
  // Print representative frames
  outfile.Printf("#Representative frames:");
  for (cluster_it C = clusters_.begin(); C != clusters_.end(); C++)
    outfile.Printf(" %i",(*C).Centroid()+1);
  outfile.Printf("\n");
  
  outfile.CloseFile();
}

// ClusterList::PrintRepFrames()
/** Print representative frame of each cluster to 1 line.
  */
void ClusterList::PrintRepFrames() {
  for (cluster_it C = clusters_.begin(); C != clusters_.end(); C++) 
    mprintf("%i ",(*C).Centroid()+1);
  mprintf("\n");
}

// -----------------------------------------------------------------------------
// ClusterList::MergeClosest()
/** Find and merge the two closest clusters.
  */
int ClusterList::MergeClosest(double epsilon, LINKAGETYPE linkage) { 
  int C1, C2;

  // Find the minimum distance between clusters. C1 will be lower than C2.
  double min = ClusterDistances_.FindMin(C1, C2);
  if (debug_>0) 
    mprintf("\tMinimum found between clusters %i and %i (%f)\n",C1,C2,min);
  // If the minimum distance is greater than epsilon we are done
  if (min > epsilon) {
    mprintf("\n\tMinimum distance (%f) is greater than epsilon (%f), clustering complete.\n",
            min, epsilon);
    return 1;
  }

  // Find the clusters in the cluster list
  // Find C1
  cluster_it C1_it = clusters_.begin();
  for (; C1_it != clusters_.end(); ++C1_it) 
  {
    if ( (*C1_it).Num() == C1 ) break;
  }
  if (C1_it == clusters_.end()) {
    mprinterr("Error: ClusterList::MergeClosest: C1 (%i) not found.\n",C1);
    return 1;
  }
  // Find C2 - start from C1 since C1 < C2
  cluster_it C2_it = C1_it;
  for (; C2_it != clusters_.end(); ++C2_it) {
    if ( (*C2_it).Num() == C2 ) break;
  }
  if (C2_it == clusters_.end()) {
    mprinterr("Error: ClusterList::MergeClosest: C2 (%i) not found.\n",C2);
    return 1;
  }

  // Merge the closest clusters, C2 -> C1
  Merge(C1_it, C2_it);
  // DEBUG
  if (debug_>1) {
    mprintf("\nAFTER MERGE of %i and %i:\n",C1,C2);
    PrintClusters();
  }
  // Remove all distances having to do with C2
  ClusterDistances_.Ignore(C2);

  // Recalculate distances between C1 and all other clusters
  switch (linkage) {
    case AVERAGELINK : calcAvgDist(C1_it); break;
    case SINGLELINK  : calcMinDist(C1_it); break;
    case COMPLETELINK: calcMaxDist(C1_it); break;
  }

  if (debug_>2) { 
    mprintf("NEW CLUSTER DISTANCES:\n");
    ClusterDistances_.PrintElements();
  }

  // Calculate the new eccentricity of C1
  //CalcEccentricity(C1_it);
  // Check if eccentricity is less than epsilon
  //if ( (*C1_it).eccentricity < epsilon) {
  //  mprintf("\tCluster %i eccentricity %lf > epsilon (%lf), clustering complete.\n",
  //          (*C1_it).num, (*C1_it).eccentricity, epsilon);
  //  return 1;
  //}

  return 0;
}

// ClusterList::Merge()
/** Merge cluster C2 into C1; remove C2. 
  */
int ClusterList::Merge(cluster_it& c1, cluster_it& c2) 
{
  // Merge C2 into C1
  (*c1).MergeFrames( *c2 );
  // Remove c2
  clusters_.erase( c2 );
  return 0;
}        

// -----------------------------------------------------------------------------
// ClusterList::calcMinDist()
/** Calculate the minimum distance between frames in cluster specified by
  * iterator C1 and frames in all other clusters.
  */
void ClusterList::calcMinDist(cluster_it& C1_it) 
{
  // All cluster distances to C1 must be recalcd.
  for (cluster_it C2_it = clusters_.begin(); 
                  C2_it != clusters_.end(); C2_it++) 
  {
    if (C2_it == C1_it) continue;
    //mprintf("\t\tRecalc distance between %i and %i:\n",C1,newc2);
    // Pick the minimum distance between newc2 and C1
    double min = DBL_MAX;
    for (ClusterNode::frame_iterator c1frames = (*C1_it).beginframe();
                                     c1frames != (*C1_it).endframe();
                                     c1frames++) 
    {
      for (ClusterNode::frame_iterator c2frames = (*C2_it).beginframe();
                                       c2frames != (*C2_it).endframe();
                                       c2frames++)
      {
        double Dist = FrameDistances_.GetElement(*c1frames, *c2frames);
        //mprintf("\t\t\tFrame %i to frame %i = %lf\n",*c1frames,*c2frames,Dist);
        if ( Dist < min ) min = Dist;
      }
    }
    //mprintf("\t\tMin distance between %i and %i: %lf\n",C1,newc2,min);
    ClusterDistances_.SetElement( (*C1_it).Num(), (*C2_it).Num(), min );
  } 
}

// ClusterList::calcMaxDist()
/** Calculate the maximum distance between frames in cluster specified by
  * iterator C1 and frames in all other clusters.
  */
void ClusterList::calcMaxDist(cluster_it& C1_it) 
{
  // All cluster distances to C1 must be recalcd.
  for (cluster_it C2_it = clusters_.begin(); 
                  C2_it != clusters_.end(); C2_it++) 
  {
    if (C2_it == C1_it) continue;
    //mprintf("\t\tRecalc distance between %i and %i:\n",C1,newc2);
    // Pick the maximum distance between newc2 and C1
    double max = -1.0;
    for (ClusterNode::frame_iterator c1frames = (*C1_it).beginframe();
                                     c1frames != (*C1_it).endframe();
                                     c1frames++) 
    {
      for (ClusterNode::frame_iterator c2frames = (*C2_it).beginframe();
                                       c2frames != (*C2_it).endframe();
                                       c2frames++)
      {
        double Dist = FrameDistances_.GetElement(*c1frames, *c2frames);
        //mprintf("\t\t\tFrame %i to frame %i = %lf\n",*c1frames,*c2frames,Dist);
        if ( Dist > max ) max = Dist;
      }
    }
    //mprintf("\t\tMax distance between %i and %i: %lf\n",C1,newc2,max);
    ClusterDistances_.SetElement( (*C1_it).Num(), (*C2_it).Num(), max );
  } 
}

// ClusterList::calcAvgDist()
/** Calculate the average distance between frames in cluster specified by
  * iterator C1 and frames in all other clusters.
  */
void ClusterList::calcAvgDist(cluster_it& C1_it) 
{
  // All cluster distances to C1 must be recalcd.
  for (cluster_it C2_it = clusters_.begin(); 
                  C2_it != clusters_.end(); C2_it++) 
  {
    if (C2_it == C1_it) continue;
    //mprintf("\t\tRecalc distance between %i and %i:\n",(*C1_it).Num(),(*C2_it).Num());
    // Pick the minimum distance between newc2 and C1
    double sumDist = 0;
    double N = 0;
    for (ClusterNode::frame_iterator c1frames = (*C1_it).beginframe();
                                     c1frames != (*C1_it).endframe();
                                     c1frames++) 
    {
      for (ClusterNode::frame_iterator c2frames = (*C2_it).beginframe();
                                       c2frames != (*C2_it).endframe();
                                       c2frames++)
      {
        double Dist = FrameDistances_.GetElement(*c1frames, *c2frames);
        //mprintf("\t\t\tFrame %i to frame %i = %lf\n",*c1frames,*c2frames,Dist);
        sumDist += Dist;
        N++;
      }
    }
    double Dist = sumDist / N;
    //mprintf("\t\tAvg distance between %i and %i: %lf\n",(*C1_it).Num(),(*C2_it).Num(),Dist);
    ClusterDistances_.SetElement( (*C1_it).Num(), (*C2_it).Num(), Dist );
  } 
}

// -----------------------------------------------------------------------------
// ClusterList::CheckEpsilon()
/** Check the eccentricity of every cluster against the given epsilon. If
  * any cluster has an eccentricity less than epsilon return true, 
  * otherwise return false.
  */
bool ClusterList::CheckEpsilon(double epsilon) {
  for (cluster_it C1_it = clusters_.begin(); C1_it != clusters_.end(); ++C1_it) 
  {
    (*C1_it).CalcEccentricity(FrameDistances_);
    if ( (*C1_it).Eccentricity() < epsilon) return true;
  }
  return false;
}

/** The Davies-Bouldin Index (DBI) is a measure of clustering merit; the 
  * smaller the DBI, the better. The DBI is defined as the average, for all 
  * clusters X, of fred, where fred(X) = max, across other clusters Y, of 
  * (Cx + Cy)/dXY ...here Cx is the average distance from points in X to the 
  * centroid, similarly Cy, and dXY is the distance between cluster centroids.
  */
double ClusterList::ComputeDBI( DataSet* dsIn, ClusterNode::RMSoptions const& rmsopt ) {
  std::vector<double> averageDist;
  averageDist.reserve( clusters_.size() );
  for (cluster_it C1 = clusters_.begin(); C1 != clusters_.end(); ++C1) {
    mprintf("AVG DISTANCES FOR CLUSTER %d:\n",(*C1).Num()); // DEBUG
    // Make sure centroid frame for this cluster is up to date
    (*C1).CalculateCentroid( dsIn, rmsopt );
    // Calculate average distance to centroid for this cluster
    averageDist.push_back( (*C1).CalcAvgToCentroid( dsIn, rmsopt ) );
    mprintf("\tCluster %i has average-distance-to-centroid %f\n", (*C1).Num(), // DEBUG
            averageDist.back()); // DEBUG
  }
  double DBITotal = 0.0;
  unsigned int nc1 = 0;
  for (cluster_it c1 = clusters_.begin(); c1 != clusters_.end(); ++c1, ++nc1) {
    double MaxFred = 0;
    unsigned int nc2 = 0;
    for (cluster_it c2 = clusters_.begin(); c2 != clusters_.end(); ++c2, ++nc2) {
      if (c1 == c2) continue;
      double Fred = averageDist[nc1] + averageDist[nc2];
      Fred /= (*c1).CentroidDist( *c2, rmsopt );
      if (Fred > MaxFred)
        MaxFred = Fred;
    }
    DBITotal += MaxFred;
  }
  return (DBITotal / (double)clusters_.size() );
}

/** The pseudo-F statistic is another measure of clustering goodness. HIGH 
  * values are GOOD. Generally, one selects a cluster-count that gives a peak 
  * in the pseudo-f statistic (or pSF, for short).
  * Formula: A/B, where A = (T - P)/(G-1), and B = P / (n-G). Here n is the 
  * number of points, G is the number of clusters, T is the total distance from
  * the all-data centroid, and P is the sum (for all clusters) of the distances
  * from the cluster centroid.
  */

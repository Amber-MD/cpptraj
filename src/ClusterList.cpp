/* Action: Clustering
 */
#include <cfloat> // DBL_MAX
#include <cmath> // sqrt
#include <vector>
#include "ClusterList.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"

// XMGRACE colors
const char ClusterList::XMGRACE_COLOR[16][12] = {
  "white", "black", "red", "green", "blue", "yellow", "brown", "grey", "violet",
  "cyan", "magenta", "orange", "indigo", "maroon", "turquoise", "darkgreen"
};

// CONSTRUCTOR
ClusterList::ClusterList() :
  debug_(0),
  maxframes_(0),
  FrameDistances_(NULL),
  Linkage_(AVERAGELINK)
{}

// ClusterList::SetDebug()
/** Set the debug level */
void ClusterList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0) mprintf("ClusterList debug set to %i\n",debug_);
}

// ClusterList::Renumber()
/** Sort clusters by size and renumber starting from 0, where cluster 0
  * is the largest.
  * NOTE: This destroys indexing into ClusterDistances.
  */
void ClusterList::Renumber() {
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
  int newNum = 0;
  for (cluster_it node = clusters_.begin();
                  node != clusters_.end(); node++) 
  {
    (*node).SetNum( newNum );
    // Find the centroid. Since FindCentroid uses FrameDistances and not
    // ClusterDistances its ok to call after sorting/renumbering.
    FindCentroid( node );
    ++newNum;
  }
}

// ClusterList::Summary()
/** Print a summary of clusters.  */
void ClusterList::Summary(const char* summaryfile) {
  CpptrajFile outfile;
  std::vector<double> distances;

  if (outfile.SetupWrite(summaryfile, debug_)) {
    mprinterr("Error: ClusterList::Summary: Could not set up file.\n");
    return;
  }
  outfile.OpenFile();

  outfile.Printf("%-8s %8s %8s %8s %8s %8s %8s\n","#Cluster","Frames","Frac",
                     "AvgDist","Stdev","Centroid","AvgCDist");
  for (cluster_it node = clusters_.begin();
                  node != clusters_.end(); node++)
  {
    // Calculate size and fraction of total size of this cluster
    int numframes = (*node).Nframes();
    float frac = (float) maxframes_;
    frac = ((float) numframes) / frac;
    // Find centroid - now done in Renumber
    //FindCentroid(node);
    // Calculate the average distance between frames in the cluster
    distances.clear();
    //int numdist = ((numframes * numframes) - numframes) / 2;
    //distances = new double[ numdist ];
    // DEBUG
    //mprintf("\tCluster %i\n",(*node).num);
    double avgdist = 0;
    //numdist = 0;
    for (ClusterNode::frame_iterator frame1 = (*node).beginframe();
                                     frame1 != (*node).endframe(); frame1++)
    {
      ClusterNode::frame_iterator frame2 = frame1;
      ++frame2;
      for (; frame2 != (*node).endframe(); frame2++) 
      {
        if (frame1==frame2) continue;
        double dist = FrameDistances_->GetElement(*frame1, *frame2);
        distances.push_back( dist );
        // DEBUG
        //mprintf("\t\tFrame %3i to %3i %8.3lf\n",*frame1,*frame2,dist);
        avgdist += dist;
        //numdist++;
      }
    }
    double sdist = 0;
    // Calculate avg and stdev of distances
    if (!distances.empty()) {
      double ndist = (double)distances.size();
      avgdist /= ndist;
      // Stdev
      for (std::vector<double>::iterator D = distances.begin();
                                         D != distances.end(); D++)
      {
        double dist = *D - avgdist;
        dist *= dist;
        sdist += dist;
      }
      sdist /= ndist;
      sdist = sqrt(sdist);
    } //else {
    //  avgdist = 0;
    //}
    // OUTPUT
    outfile.Printf("%8i %8i %8.3f %8.3lf %8.3lf %8i %8.3lf\n",
                   (*node).Num(), numframes, frac, avgdist, sdist,
                   (*node).Centroid()+1, (*node).AvgDist() );
    //delete[] distances;
  } // END loop over clusters

  outfile.CloseFile();
}

// ClusterList::Summary_Half
/** Print a summary of the first half of the data to the second half.
  */
void ClusterList::Summary_Half(const char* summaryfile) {
  CpptrajFile outfile;
/*  int numInFirstHalf, numInSecondHalf;
  int numframes, half, color;
  float frac, frac1, frac2;*/

  if (outfile.SetupWrite(summaryfile, debug_)) {
    mprinterr("Error: ClusterList::Summary_Half: Could not set up file.\n");
    return;
  }
  outfile.OpenFile();

  // Calculate halfway point
  int half = maxframes_ / 2;
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
    float frac = (float) maxframes_;
    frac = ((float) numframes) / frac;
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

// ClusterList::AddCluster()
/** Add a cluster made up of frames specified by the given framelist to 
  * the list.
  */
int ClusterList::AddCluster( std::list<int>& framelistIn, int numIn ) {
  // Update number of frames
  // NOTE: Now done in initialization
  //maxframes_ += framelistIn.size();
  clusters_.push_back( ClusterNode( framelistIn, numIn ) );
  return 0;
}

// ClusterList::Initialize()
/** Given a triangle matrix containing the distances between all frames,
  * set up the initial distances between clusters.
  * Should be called before any clustering is performed. 
  */
void ClusterList::Initialize(TriangleMatrix *matrixIn) {
  FrameDistances_ = matrixIn;
  maxframes_ = matrixIn->Nrows();
  ClusterDistances_.Setup( clusters_.size() );
  // Build initial cluster distances
  if (Linkage_==AVERAGELINK) {
    for (cluster_it C1_it = clusters_.begin(); 
                    C1_it != clusters_.end(); C1_it++) 
      calcAvgDist(C1_it);
  } else if (Linkage_==SINGLELINK) {
    for (cluster_it C1_it = clusters_.begin();
                    C1_it != clusters_.end(); C1_it++) 
      calcMinDist(C1_it);
  } else if (Linkage_==COMPLETELINK) {
    for (cluster_it C1_it = clusters_.begin(); 
                    C1_it != clusters_.end(); C1_it++) 
      calcMaxDist(C1_it);
  }
  if (debug_ > 1) {
    mprintf("CLUSTER: INITIAL CLUSTER DISTANCES:\n");
    ClusterDistances_.PrintElements();
  }
}
 
// ClusterList::PrintClusters()
/** Print list of clusters and frame numbers belonging to each cluster.
  */
void ClusterList::PrintClusters() {
  mprintf("CLUSTER: %u clusters, %i frames.\n", clusters_.size(),maxframes_);
  for (cluster_it C = clusters_.begin(); C != clusters_.end(); C++) {
    mprintf("\t%8i : ",(*C).Num());
    for (ClusterNode::frame_iterator fnum = (*C).beginframe();
                                     fnum != (*C).endframe();
                                     fnum++)
    {
      mprintf("%i,",(*fnum)+1);
    }
    mprintf("\n");
  }
}

// ClusterList::PrintClustersToFile()
/** Print list of clusters in a style similar to ptraj; each cluster is
  * given a line maxframes characters long, with X for each frame that is
  * in the clusters and . for all other frames. Also print out the
  * representative frame numbers.
  */
void ClusterList::PrintClustersToFile(const char* filename) {
  CpptrajFile outfile;
  std::string buffer;
  
  if ( outfile.SetupWrite(filename,debug_) ) {
    mprinterr("Error: PrintClustersToFile: Could not set up file %s\n",
              filename);
    return;
  }
  outfile.OpenFile();
  outfile.Printf("#Clustering: %u clusters %i frames\n",
                 clusters_.size(), maxframes_);
  for (cluster_it C1_it = clusters_.begin(); 
                  C1_it != clusters_.end(); C1_it++)
  {
    buffer.clear();
    buffer.resize(maxframes_, '.');
    for (ClusterNode::frame_iterator frame1 = (*C1_it).beginframe();
                                     frame1 != (*C1_it).endframe();
                                     frame1++)
    {
      buffer[ *frame1 ] = 'X';
    }
    buffer += '\n';
    outfile.Write((char*)buffer.c_str(), buffer.size());
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
int ClusterList::MergeClosest(double epsilon) { 
  int C1, C2;

  // Find the minimum distance between clusters. C1 will be lower than C2.
  double min = ClusterDistances_.FindMin(&C1, &C2);
  if (debug_>0) 
    mprintf("\tMinimum found between clusters %i and %i (%lf)\n",C1,C2,min);
  // If the minimum distance is greater than epsilon we are done
  if (min > epsilon) {
    mprintf("\tMinimum distance is greater than epsilon (%lf), clustering complete.\n", epsilon);
    return 1;
  }

  // Find the clusters in the cluster list
  // Find C1
  cluster_it C1_it = clusters_.begin();
  for (; C1_it != clusters_.end(); C1_it++) 
  {
    if ( (*C1_it).Num() == C1 ) break;
  }
  if (C1_it == clusters_.end()) {
    mprinterr("Error: ClusterList::MergeClosest: C1 (%i) not found.\n",C1);
    return 1;
  }
  // Find C2 - start from C1 since C1 < C2
  cluster_it C2_it = C1_it;
  for (; C2_it != clusters_.end(); C2_it++) {
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
  if (Linkage_==AVERAGELINK)
    calcAvgDist(C1_it);
  else if (Linkage_==SINGLELINK)
    calcMinDist(C1_it);
  else if (Linkage_==COMPLETELINK)
    calcMaxDist(C1_it);

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
        double Dist = FrameDistances_->GetElement(*c1frames, *c2frames);
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
        double Dist = FrameDistances_->GetElement(*c1frames, *c2frames);
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
        double Dist = FrameDistances_->GetElement(*c1frames, *c2frames);
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

// ClusterList::FindCentroid()
/** Find the frame in the given cluster that is the centroid, i.e. has the
  * lowest cumulative distance to every other point in the cluster.
  */
void ClusterList::FindCentroid(cluster_it& C1_it) {
  double mindist = DBL_MAX;
  int minframe = -1;
  for (ClusterNode::frame_iterator frame1 = (*C1_it).beginframe(); 
                                   frame1 != (*C1_it).endframe();
                                   frame1++)
  {
    double cdist = 0;
    for (ClusterNode::frame_iterator frame2 = (*C1_it).beginframe();
                                     frame2 != (*C1_it).endframe();
                                     frame2++)
    {
      if (frame1==frame2) continue;
      cdist += FrameDistances_->GetElement(*frame1, *frame2);
    }
    if (cdist < mindist) {
      mindist = cdist;
      minframe = (*frame1);
    }
  }
  if (minframe==-1) {
    mprinterr("Error: FindCentroid: Cluster %i could not determine centroid frame.\n", 
              (*C1_it).Num());
    return;
  }
  (*C1_it).SetCentroid( minframe );
}

// ClusterList::CalcEccentricity()
/** Calculate the eccentricity of the given cluster, i.e. the largest distance
  * between any two points in the cluster.
  */
void ClusterList::CalcEccentricity(cluster_it& C1_it) {
  double maxdist = 0;
  ClusterNode::frame_iterator frame1_end = (*C1_it).endframe();
  --frame1_end;

  for (ClusterNode::frame_iterator frame1 = (*C1_it).beginframe();
                                   frame1 != frame1_end;
                                   frame1++)
  {
    ClusterNode::frame_iterator frame2 = frame1;
    ++frame2;
    for (; frame2 != (*C1_it).endframe(); frame2++) {
      double fdist = FrameDistances_->GetElement(*frame1, *frame2);
      if (fdist > maxdist) maxdist = fdist;
    }
  } 
  (*C1_it).SetEccentricity( maxdist );
}

// ClusterList::CheckEpsilon()
/** Check the eccentricity of every cluster against the given epsilon. If
  * any cluster has an eccentricity less than epsilon return true, 
  * otherwise return false.
  */
bool ClusterList::CheckEpsilon(double epsilon) {
  for (cluster_it C1_it = clusters_.begin();
                  C1_it != clusters_.end(); C1_it++) 
  {
    CalcEccentricity( C1_it );
    if ( (*C1_it).Eccentricity() < epsilon) return true;
  }
  return false;
}


/* Action: Clustering
 */
#include <cfloat>
#include <cmath>
#include <cstring> // memset
#include "ClusterList.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"

// XMGRACE colors
// NOTE: Should this be somewhere else?
const char XMGRACE_COLOR[16][12] = {
"white",
"black",
"red",
"green",
"blue",
"yellow",
"brown",
"grey",
"violet",
"cyan",
"magenta",
"orange",
"indigo",
"maroon",
"turquoise",
"darkgreen" };

// CONSTRUCTOR
ClusterList::ClusterList() {
  debug=0;
  maxframes=0;
  Linkage = AVERAGELINK;
}

// ClusterList::SetDebug()
/** Set the debug level */
void ClusterList::SetDebug(int debugIn) {
  debug = debugIn;
  if (debug>0) mprintf("ClusterList debug set to %i\n",debug);
}

// ClusterList::Renumber()
/** Sort clusters by size and renumber starting from 0, where cluster 0
  * is the largest.
  * NOTE: This destroys indexing into ClusterDistances.
  */
void ClusterList::Renumber() {
  int newNum = 0;
  double numdist;

  // Before clusters are renumbered, calculate the average distance of 
  // this cluster to every other cluster
  numdist = (double) (clusters.size() - 1);
  for (std::list<clusterNode>::iterator node = clusters.begin();
                                        node != clusters.end();
                                        node++)
  {
    (*node).avgclusterdist = 0;
    for (std::list<clusterNode>::iterator node2 = clusters.begin();
                                          node2 != clusters.end();
                                          node2++)
    {
      if (node == node2) continue;
      //mprintf("DBG:\t\t%i to %i %lf\n",(*node).num, (*node2).num, 
      //        ClusterDistances.GetElement( (*node).num, (*node2).num ));
      (*node).avgclusterdist += ClusterDistances.GetElement( (*node).num, (*node2).num );
    }
    (*node).avgclusterdist /= numdist;
    //mprintf("DBG:\tCluster %i avg dist = %lf\n",(*node).num,avgclusterdist);
  }
  
  clusters.sort( cluster_cmp() );
  for (std::list<clusterNode>::iterator node = clusters.begin();
                                        node != clusters.end();
                                        node++) 
  {
    (*node).num = newNum;
    // Sort the frame lists for good measure
    (*node).frameList.sort();
    // Find the centroid. Since FindCentroid uses FrameDistances and not
    // ClusterDistances its ok to call after sorting/renumbering.
    FindCentroid( node );
    newNum++;
  }
}

// ClusterList::Summary()
/** Print a summary of clusters.  */
void ClusterList::Summary(char *summaryfile) {
  CpptrajFile outfile;
  int numframes,numdist;
  float frac;
  double dist,avgdist,sdist,*distances;

  if (outfile.SetupFile(summaryfile, WRITE, 0)) {
    mprinterr("Error: ClusterList::Summary: Could not set up file.\n");
    return;
  }
  outfile.OpenFile();

  outfile.IO->Printf("%-8s %8s %8s %8s %8s %8s %8s\n","#Cluster","Frames","Frac",
                     "AvgDist","Stdev","Centroid","AvgCDist");
  for (std::list<clusterNode>::iterator node = clusters.begin();
                                        node != clusters.end();
                                        node++)
  {
    // Calculate size and fraction of total size of this cluster
    numframes = (*node).frameList.size();
    frac = (float) maxframes;
    frac = ((float) numframes) / frac;
    // Find centroid - now done in Renumber
    //FindCentroid(node);
    // Calculate the average distance between frames in the cluster
    numdist = ((numframes * numframes) - numframes) / 2;
    distances = new double[ numdist ];
    // DEBUG
    //mprintf("\tCluster %i\n",(*node).num);
    avgdist = 0;
    numdist = 0;
    for (std::list<int>::iterator frame1 = (*node).frameList.begin();
                                  frame1 != (*node).frameList.end();
                                  frame1++)
    {
      std::list<int>::iterator frame2 = frame1;
      frame2++;
      for (; frame2 != (*node).frameList.end(); frame2++) 
      {
        if (frame1==frame2) continue;
        dist = FrameDistances->GetElement(*frame1,*frame2);
        distances[numdist] = dist;
        // DEBUG
        //mprintf("\t\tFrame %3i to %3i %8.3lf\n",*frame1,*frame2,dist);
        avgdist += dist;
        numdist++;
      }
    }
    if (numdist > 0) {
      avgdist /= ((double) numdist);
      // Stdev
      sdist = 0;
      for (int N=0; N < numdist; N++) {
        dist = distances[N] - avgdist;
        dist *= dist;
        sdist += dist;
      }
      sdist /= ((double) numdist);
      sdist = sqrt(sdist);
    } else {
      avgdist = 0;
      sdist = 0;
    }
    // OUTPUT
    outfile.IO->Printf("%8i %8i %8.3f %8.3lf %8.3lf %8i %8.3lf\n",(*node).num,numframes,
                       frac,avgdist,sdist,(*node).centroid+1,(*node).avgclusterdist);
    delete[] distances;
  }

  outfile.CloseFile();
}

// ClusterList::Summary_Half
/** Print a summary of the first half of the data to the second half.
  */
void ClusterList::Summary_Half(char *summaryfile) {
  CpptrajFile outfile;
  int numInFirstHalf, numInSecondHalf;
  int numframes, half, color;
  float frac, frac1, frac2;

  if (outfile.SetupFile(summaryfile, WRITE, debug)) {
    mprinterr("Error: ClusterList::Summary_Half: Could not set up file.\n");
    return;
  }
  outfile.OpenFile();

  // Calculate halfway point
  half = maxframes / 2;
  // xmgrace color
  color=1;

  outfile.IO->Printf("#%-7s %8s %6s %2s %10s %8s %8s %6s %6s\n", "Cluster", "Total",
                     "Frac", "C#", "Color", "NumIn1st", "NumIn2nd","Frac1","Frac2");
  for (std::list<clusterNode>::iterator node = clusters.begin();
                                        node != clusters.end();
                                        node++)
  {
    // Calculate size and fraction of total size of this cluster
    numframes = (*node).frameList.size();
    frac = (float) maxframes;
    frac = ((float) numframes) / frac;
    numInFirstHalf=0;
    numInSecondHalf=0;
    // DEBUG
    //mprintf("\tCluster %i\n",(*node).num);
    // Count how many frames are in the first half and how many 
    // are in the second half.
    for (std::list<int>::iterator frame1 = (*node).frameList.begin();
                                  frame1 != (*node).frameList.end();
                                  frame1++)
    {
      if (*frame1 < half)
        numInFirstHalf++;
      else
        numInSecondHalf++;
    }
    frac1 = (float) numframes;
    frac1 = ((float) numInFirstHalf) / frac1;
    frac2 = (float) numframes;
    frac2 = ((float) numInSecondHalf) / frac2;
    outfile.IO->Printf("%-8i %8i %6.2f %2i %10s %8i %8i %6.2f %6.2f\n",(*node).num,numframes,
                       frac,color,XMGRACE_COLOR[color],numInFirstHalf,numInSecondHalf,
                       frac1,frac2);
    if (color<15) color++;
  }
  outfile.CloseFile();
}

// ClusterList::AddCluster()
/** Add a cluster made up of frames specified by the given framelist to 
  * the list.
  */
int ClusterList::AddCluster( std::list<int> *framelistIn, int numIn  ) {
  clusterNode CN;

  CN.frameList = *framelistIn;
  CN.num = numIn;
  // Set initial centroid to front, even though that will probably be wrong
  // when number of frames in the list > 1
  CN.centroid = framelistIn->front();
  maxframes += CN.frameList.size();

  clusters.push_back(CN);

  return 0;
}

// ClusterList::Initialize()
/** Given a triangle matrix containing the distances between all frames,
  * set up the initial distances between clusters.
  * Should be called before any clustering is performed. 
  */
void ClusterList::Initialize(TriangleMatrix *matrixIn) {
  std::list<clusterNode>::iterator C1_it;

  FrameDistances = matrixIn;
  ClusterDistances.Setup( FrameDistances->Nrows() );
  // Build initial cluster distances
  if (Linkage==AVERAGELINK) {
    for (C1_it = clusters.begin(); C1_it != clusters.end(); C1_it++) 
      calcAvgDist(C1_it);
  } else if (Linkage==SINGLELINK) {
    for (C1_it = clusters.begin(); C1_it != clusters.end(); C1_it++) 
      calcMinDist(C1_it);
  } else if (Linkage==COMPLETELINK) {
    for (C1_it = clusters.begin(); C1_it != clusters.end(); C1_it++) 
      calcMaxDist(C1_it);
  }
}
    
// ClusterList::PrintClusters()
/** Print list of clusters and frame numbers belonging to each cluster.
  */
void ClusterList::PrintClusters() {
  mprintf("CLUSTER: %u clusters, %i frames.\n", clusters.size(),maxframes);
  for (std::list<clusterNode>::iterator C = clusters.begin(); C != clusters.end(); C++) {
    mprintf("\t%8i : ",(*C).num);
    for (std::list<int>::iterator fnum = (*C).frameList.begin();
                                  fnum != (*C).frameList.end();
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
void ClusterList::PrintClustersToFile(char *filename) {
  CpptrajFile outfile;
  char *buffer;
  int buffersize = maxframes + 1; // + 1 for newline char
  
  buffer = new char[ buffersize ]; 
  buffer[maxframes]='\n';
  if ( outfile.SetupFile(filename,WRITE,debug) ) {
    mprinterr("Error: ClusterList::PrintClustersToFile: Could not set up file %s\n",filename);
    return;
  }
  outfile.OpenFile();
  outfile.IO->Printf("#Clustering: %u clusters %i frames\n",clusters.size(),maxframes);
  for (std::list<clusterNode>::iterator C1_it = clusters.begin(); 
       C1_it != clusters.end();
       C1_it++)
  {
    memset(buffer,'.',maxframes);
    for (std::list<int>::iterator frame1 = (*C1_it).frameList.begin();
                                  frame1 != (*C1_it).frameList.end();
                                  frame1++)
    {
      buffer[ *frame1 ]='X';
    }
    outfile.IO->Write(buffer, sizeof(char), buffersize);
  }
  // Print representative frames
  outfile.IO->Printf("#Representative frames:");
  for (std::list<clusterNode>::iterator C = clusters.begin(); C != clusters.end(); C++)
    outfile.IO->Printf(" %i",(*C).centroid+1);
  outfile.IO->Printf("\n");
  
  outfile.CloseFile();
  delete[] buffer;
}

// ClusterList::PrintRepFrames()
/** Print representative frame of each cluster to 1 line.
  */
void ClusterList::PrintRepFrames() {
  for (std::list<clusterNode>::iterator C = clusters.begin(); C != clusters.end(); C++) 
    mprintf("%i ",(*C).centroid+1);
  mprintf("\n");
}

// ClusterList::GetClusterIt
/** Return an iterator to the specified cluster.
  */
std::list<ClusterList::clusterNode>::iterator ClusterList::GetClusterIt(int C1) {
  std::list<clusterNode>::iterator c1;
  // Find C1
  for (c1 = clusters.begin(); c1 != clusters.end(); c1++) {
    if ( (*c1).num == C1 ) break;
  }
  if (c1 == clusters.end()) {
    mprinterr("Error: ClusterList::Merge: C1 (%i) not found.\n",C1);
    //return 1;
  }
  return c1;
}

// ClusterList::MergeClosest()
/** Find and merge the two closest clusters.
  */
int ClusterList::MergeClosest(double epsilon) { 
  double min;
  int C1, C2;
  std::list<clusterNode>::iterator C1_it;
  std::list<clusterNode>::iterator C2_it;

  // Find the minimum distance between clusters. C1 will be lower than C2.
  min = ClusterDistances.FindMin(&C1, &C2);
  if (debug>0) mprintf("\tMinimum found between clusters %i and %i (%lf)\n",C1,C2,min);
  // If the minimum distance is greater than epsilon we are done
  if (min > epsilon) {
    mprintf("\tMinimum distance is greater than epsilon (%lf), clustering complete.\n",
            epsilon);
    return 1;
  }

  // Find the clusters in the cluster list
  // Find C1
  for (C1_it = clusters.begin(); C1_it != clusters.end(); C1_it++) {
    if ( (*C1_it).num == C1 ) break;
  }
  if (C1_it == clusters.end()) {
    mprinterr("Error: ClusterList::MergeClosest: C1 (%i) not found.\n",C1);
    return 1;
  }
  // Find C2 - start from C1 since C1 < C2
  for (C2_it = C1_it; C2_it != clusters.end(); C2_it++) {
    if ( (*C2_it).num == C2 ) break;
  }
  if (C2_it == clusters.end()) {
    mprinterr("Error: ClusterList::MergeClosest: C2 (%i) not found.\n",C2);
    return 1;
  }

  // Merge the closest clusters, C2 -> C1
  Merge(C1_it,C2_it);
  // DEBUG
  if (debug>1) {
    mprintf("\nAFTER MERGE of %i and %i:\n",C1,C2);
    PrintClusters();
  }
  // Remove all distances having to do with C2
  ClusterDistances.Ignore(C2);

  // Recalculate distances between C1 and all other clusters
  if (Linkage==AVERAGELINK)
    calcAvgDist(C1_it);
  else if (Linkage==SINGLELINK)
    calcMinDist(C1_it);
  else if (Linkage==COMPLETELINK)
    calcMaxDist(C1_it);

  if (debug>2) { 
    mprintf("NEW CLUSTER DISTANCES:\n");
    ClusterDistances.PrintElements();
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
int ClusterList::Merge(std::list<ClusterList::clusterNode>::iterator c1, 
                       std::list<ClusterList::clusterNode>::iterator c2) 
{
  // Merge C2 into C1
  std::list<int>::iterator frameit = (*c1).frameList.begin();
  (*c1).frameList.splice( frameit, (*c2).frameList );
  // Remove c2
  clusters.erase( c2 );

  return 0;
}        

// ClusterList::calcMinDist()
/** Calculate the minimum distance between frames in cluster specified by
  * iterator C1 and frames in all other clusters.
  */
void ClusterList::calcMinDist(std::list<ClusterList::clusterNode>::iterator C1_it) 
{
  double min, Dist;
  std::list<clusterNode>::iterator C2_it;

  // All cluster distances to C1 must be recalcd.
  for (C2_it = clusters.begin(); C2_it != clusters.end(); C2_it++) {
    if (C2_it == C1_it) continue;
    //mprintf("\t\tRecalc distance between %i and %i:\n",C1,newc2);
    // Pick the minimum distance between newc2 and C1
    min = DBL_MAX;
    for (std::list<int>::iterator c1frames = (*C1_it).frameList.begin();
                                  c1frames != (*C1_it).frameList.end();
                                  c1frames++) 
    {
      for (std::list<int>::iterator c2frames = (*C2_it).frameList.begin();
                                    c2frames != (*C2_it).frameList.end();
                                    c2frames++)
      {
        Dist = FrameDistances->GetElement(*c1frames, *c2frames);
        //mprintf("\t\t\tFrame %i to frame %i = %lf\n",*c1frames,*c2frames,Dist);
        if ( Dist < min ) min = Dist;
      }
    }
    //mprintf("\t\tMin distance between %i and %i: %lf\n",C1,newc2,min);
    ClusterDistances.SetElement( (*C1_it).num, (*C2_it).num, min );
  } 
}

// ClusterList::calcMaxDist()
/** Calculate the maximum distance between frames in cluster specified by
  * iterator C1 and frames in all other clusters.
  */
void ClusterList::calcMaxDist(std::list<ClusterList::clusterNode>::iterator C1_it) 
{
  double max, Dist;
  std::list<clusterNode>::iterator C2_it;

  // All cluster distances to C1 must be recalcd.
  for (C2_it = clusters.begin(); C2_it != clusters.end(); C2_it++) {
    if (C2_it == C1_it) continue;
    //mprintf("\t\tRecalc distance between %i and %i:\n",C1,newc2);
    // Pick the maximum distance between newc2 and C1
    max = -1.0;
    for (std::list<int>::iterator c1frames = (*C1_it).frameList.begin();
                                  c1frames != (*C1_it).frameList.end();
                                  c1frames++) 
    {
      for (std::list<int>::iterator c2frames = (*C2_it).frameList.begin();
                                    c2frames != (*C2_it).frameList.end();
                                    c2frames++)
      {
        Dist = FrameDistances->GetElement(*c1frames, *c2frames);
        //mprintf("\t\t\tFrame %i to frame %i = %lf\n",*c1frames,*c2frames,Dist);
        if ( Dist > max ) max = Dist;
      }
    }
    //mprintf("\t\tMax distance between %i and %i: %lf\n",C1,newc2,max);
    ClusterDistances.SetElement( (*C1_it).num, (*C2_it).num, max );
  } 
}

// ClusterList::calcAvgDist()
/** Calculate the average distance between frames in cluster specified by
  * iterator C1 and frames in all other clusters.
  */
void ClusterList::calcAvgDist(std::list<ClusterList::clusterNode>::iterator C1_it) 
{
  double N, Dist, sumDist;
  std::list<clusterNode>::iterator C2_it;

  // All cluster distances to C1 must be recalcd.
  for (C2_it = clusters.begin(); C2_it != clusters.end(); C2_it++) {
    if (C2_it == C1_it) continue;
    //mprintf("\t\tRecalc distance between %i and %i:\n",(*C1_it).num,(*C2_it).num);
    // Pick the minimum distance between newc2 and C1
    sumDist=0.0;
    N=0.0;
    for (std::list<int>::iterator c1frames = (*C1_it).frameList.begin();
                                  c1frames != (*C1_it).frameList.end();
                                  c1frames++) 
    {
      for (std::list<int>::iterator c2frames = (*C2_it).frameList.begin();
                                    c2frames != (*C2_it).frameList.end();
                                    c2frames++)
      {
        Dist = FrameDistances->GetElement(*c1frames, *c2frames);
        //mprintf("\t\t\tFrame %i to frame %i = %lf\n",*c1frames,*c2frames,Dist);
        sumDist += Dist;
        N++;
      }
    }
    Dist = sumDist / N;
    //mprintf("\t\tAvg distance between %i and %i: %lf\n",(*C1_it).num,(*C2_it).num,Dist);
    ClusterDistances.SetElement( (*C1_it).num, (*C2_it).num, Dist );
  } 
}

// ClusterList::FindCentroid()
/** Find the frame in the given cluster that is the centroid, i.e. has the
  * lowest cumulative distance to every other point in the cluster.
  */
void ClusterList::FindCentroid(std::list<ClusterList::clusterNode>::iterator C1_it) {
  double mindist = DBL_MAX;
  double cdist;
  int minframe = -1;
  for (std::list<int>::iterator frame1 = (*C1_it).frameList.begin(); 
                                frame1 != (*C1_it).frameList.end();
                                frame1++)
  {
    cdist = 0;
    for (std::list<int>::iterator frame2 = (*C1_it).frameList.begin();
                                  frame2 != (*C1_it).frameList.end();
                                  frame2++)
    {
      if (frame1==frame2) continue;
      cdist += FrameDistances->GetElement(*frame1, *frame2);
    }
    if (cdist < mindist) {
      mindist = cdist;
      minframe = (*frame1);
    }
  }
  if (minframe==-1) {
    mprinterr("Error: ClusterList::FindCentroid: Cluster %i could not determine centroid frame.\n",
              (*C1_it).num);
    return;
  }
  (*C1_it).centroid = minframe;
}

// ClusterList::CalcEccentricity()
/** Calculate the eccentricity of the given cluster, i.e. the largest distance
  * between any two points in the cluster.
  */
void ClusterList::CalcEccentricity(std::list<ClusterList::clusterNode>::iterator C1_it) {
  double maxdist = 0;
  double fdist;
  std::list<int>::iterator frame1_end = (*C1_it).frameList.end();
  frame1_end--;

  for (std::list<int>::iterator frame1 = (*C1_it).frameList.begin();
                                frame1 != frame1_end;
                                frame1++)
  {
    std::list<int>::iterator frame2 = frame1;
    frame2++;
    for (; frame2 != (*C1_it).frameList.end(); frame2++) {
      fdist = FrameDistances->GetElement(*frame1, *frame2);
      if (fdist > maxdist) maxdist = fdist;
    }
  } 
  (*C1_it).eccentricity = maxdist;
}

// ClusterList::CheckEpsilon()
/** Check the eccentricity of every cluster against the given epsilon. If
  * any cluster has an eccentricity less than epsilon return true, 
  * otherwise return false.
  */
bool ClusterList::CheckEpsilon(double epsilon) {
  for (std::list<clusterNode>::iterator C1_it = clusters.begin();
       C1_it != clusters.end();
       C1_it++) 
  {
    CalcEccentricity( C1_it );
    if ( (*C1_it).eccentricity < epsilon) return true;
  }
  return false;
}

// ClusterList::Begin()
/** Place current cluster at beginning of list.
  */
void ClusterList::Begin() {
  currentCluster = clusters.begin();
}

// ClusterList::End()
/** Return true if current cluster is at the end of list.
  */
bool ClusterList::End() {
  if (currentCluster == clusters.end()) return true;
  return false;
}

// ClusterList::NextCluster()
/** Advance current cluster to the next cluster.
  */
void ClusterList::NextCluster() {
  currentCluster++;
}

// ClusterList::CurrentNum()
/** Return number of the current cluster.  */
int ClusterList::CurrentNum() {
  return (*currentCluster).num;
}

// ClusterList::CurrentCentroid()
/** Return frame number of centroid of current cluster.  */
int ClusterList::CurrentCentroid() {
  // FindCentroid is now called in Renumber
  //FindCentroid(currentCluster);
  return (*currentCluster).centroid;
}

// ClusterList::CurrentFrameBegin()
/** Return iterator to the beginning of the current clusters framelist.
  */
std::list<int>::iterator ClusterList::CurrentFrameBegin() {
  return (*currentCluster).frameList.begin();
}

// ClusterList::CurrentFrameEnd()
/** Return iterator to the end of the current clusters framelist.
  */
std::list<int>::iterator ClusterList::CurrentFrameEnd() {
  return (*currentCluster).frameList.end();
}


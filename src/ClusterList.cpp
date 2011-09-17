/* Action: Clustering
 */
#include "ClusterList.h"
#include "CpptrajStdio.h"
#include "PtrajFile.h"
#include <cfloat>

// CONSTRUCTOR
ClusterList::ClusterList() {
  debug=0;
  maxframes=0;
  Linkage = AVERAGELINK;
}

// DESTRUCTOR
ClusterList::~ClusterList() {
}

/* ClusterList::Renumber()
 * Sort clusters by size and renumber starting from 0, where cluster 0
 * is the largest.
 */
void ClusterList::Renumber() {
  int newNum = 0;

  clusters.sort( cluster_cmp() );
  for (std::list<clusterNode>::iterator node = clusters.begin();
                                        node != clusters.end();
                                        node++) 
  {
    (*node).num = newNum;
    // Sort the frame lists for good measure
    (*node).frameList.sort();
    newNum++;
  }
}

/* ClusterList::Summary()
 * Print a summary of clusters.
 */
void ClusterList::Summary(char *summaryfile) {
  PtrajFile outfile;
  int numframes;
  float frac;

  if (outfile.SetupFile(summaryfile, WRITE, 0)) {
    mprinterr("Error: ClusterList::Summary: Could not set up file.\n");
    return;
  }
  outfile.OpenFile();

  outfile.IO->Printf("%-8s %8s %8s\n","#Cluster","Frames","Frac");
  for (std::list<clusterNode>::iterator node = clusters.begin();
                                        node != clusters.end();
                                        node++)
  {
    numframes = (*node).frameList.size();
    frac = (float) maxframes;
    frac = ((float) numframes) / frac;
    outfile.IO->Printf("%8i %8i %8.3f\n",(*node).num,numframes,frac);
  }

  outfile.CloseFile();
}

/* ClusterList::AddCluster()
 * Add a cluster made up of frames specified by the given framelist to 
 * the list.
 */
int ClusterList::AddCluster( std::list<int> *framelistIn, int numIn  ) {
  clusterNode CN;

  CN.frameList = *framelistIn;
  CN.num = numIn;
  maxframes += CN.frameList.size();

  clusters.push_back(CN);

  return 0;
}

/* ClusterList::Initialize()
 * Given a triangle matrix containing the distances between all frames,
 * set up the initial distances between clusters.
 * Should be called before any clustering is performed. 
 */
void ClusterList::Initialize(TriangleMatrix *matrixIn) {
  std::list<clusterNode>::iterator C1_it;
  std::list<clusterNode>::iterator C2_it;

  FrameDistances = matrixIn;
  ClusterDistances.Setup( FrameDistances->Nrows() );
  // Build initial cluster distances
  // NOTE: Modify the calcXDist routines to take an iterator.
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
    
/* ClusterList::PrintClusters()
 * Print list of clusters and frame numbers belonging to each cluster.
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

/* ClusterList::GetClusterIt
 * Return an iterator to the specified cluster.
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

/* ClusterList::MergeClosest()
 * Find and merge the two closest clusters.
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

  // Merge the closest clusters
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
  return 0;
}

/* ClusterList::Merge()
 * Merge cluster C2 into C1; remove C2. 
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

/* ClusterList::calcMinDist()
 * Calculate the minimum distance between frames in cluster specified by
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

/* ClusterList::calcMaxDist()
 * Calculate the maximum distance between frames in cluster specified by
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

/* ClusterList::calcAvgDist()
 * Calculate the average distance between frames in cluster specified by
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

/* ClusterList::Cnumvtime()
 * Put cluster number v frame into dataset.
 */
void ClusterList::Cnumvtime(DataSet *cnumvtimeIn) {
  int cnum;

  for (std::list<clusterNode>::iterator node = clusters.begin();
                                        node != clusters.end();
                                        node++)
  {
    for (std::list<int>::iterator frame = (*node).frameList.begin();
                                  frame != (*node).frameList.end();
                                  frame++) {
      cnum = (*node).num;
      cnumvtimeIn->Add( *frame, &cnum );
    }
  }
}


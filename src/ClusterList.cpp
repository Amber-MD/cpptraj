/* Action: Clustering
 */
#include "ClusterList.h"
#include "CpptrajStdio.h"
#include "PtrajFile.h"
#include <cfloat>

// CONSTRUCTOR
ClusterList::ClusterList() {
  maxframes=0;
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

/* ClusterList::GetCluster
 * Return an iterator to the specified cluster.
 */
std::list<ClusterList::clusterNode>::iterator ClusterList::GetCluster(int C1) {
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
void ClusterList::calcMinDist(std::list<ClusterList::clusterNode>::iterator C1_it,
                              TriangleMatrix *FrameDistances,
                              TriangleMatrix *ClusterDistances) 
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
    //mprintf("\t\tNew distance between %i and %i: %lf\n",C1,newc2,min);
    ClusterDistances->SetElement( (*C1_it).num, (*C2_it).num, min );
  } 
}

/* ClusterList::calcAvgDist()
 * Calculate the average distance between frames in cluster specified by
 * iterator C1 and frames in all other clusters.
 */
void ClusterList::calcAvgDist(std::list<ClusterList::clusterNode>::iterator C1_it,
                              TriangleMatrix *FrameDistances,
                              TriangleMatrix *ClusterDistances) 
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
    ClusterDistances->SetElement( (*C1_it).num, (*C2_it).num, Dist );
  } 
}

/* ClusterList::Cnumvtime()
 * Put cluster number v frame into dataset.
 */
void ClusterList::Cnumvtime(DataSet *cnumvtimeIn) {
  // Temp array
  int *CVT = new int[ maxframes ];

  for (std::list<clusterNode>::iterator node = clusters.begin();
                                        node != clusters.end();
                                        node++)
  {
    for (std::list<int>::iterator frame = (*node).frameList.begin();
                                  frame != (*node).frameList.end();
                                  frame++) {
      CVT[ *frame ] = (*node).num;
    }
  }

  for (int N=0; N < maxframes; N++) 
    cnumvtimeIn->Add(N,CVT + N);
  delete[] CVT;
}



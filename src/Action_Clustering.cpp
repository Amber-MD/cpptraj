/* Action: Clustering
 */
#include "Action_Clustering.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#include <cfloat>
// ============================================================================
/// Class: ClusterList
#include <list>
/// Store information on clusters
class ClusterList {
  public :
    // Store individual cluster info; frame numbers, centroid, etc.
    struct clusterNode {
      std::list<int> frameList; // List of frames in the cluster
      int num;                  // Cluster number (index in Distances)
    };
    std::list<clusterNode> clusters;

    ClusterList();
    ~ClusterList();
   
    int AddCluster(std::list<int> *, int);
    void PrintClusters();
    std::list<clusterNode>::iterator GetCluster(int);
    int Merge(std::list<clusterNode>::iterator,std::list<clusterNode>::iterator);
    void calcMinDist(std::list<ClusterList::clusterNode>::iterator,
                     TriangleMatrix *, TriangleMatrix *);
    void calcAvgDist(std::list<ClusterList::clusterNode>::iterator,
                     TriangleMatrix *, TriangleMatrix *);
};

// CONSTRUCTOR
ClusterList::ClusterList() {
}

// DESTRUCTOR
ClusterList::~ClusterList() {
}

/* ClusterList::AddCluster()
 * Add a cluster made up of frames specified by the given framelist to 
 * the list.
 */
int ClusterList::AddCluster( std::list<int> *framelistIn, int numIn  ) {
  clusterNode CN;

  CN.frameList = *framelistIn;
  CN.num = numIn;

  clusters.push_back(CN);

  return 0;
}

/* ClusterList::PrintClusters()
 * Print list of clusters and frame numbers belonging to each cluster.
 */
void ClusterList::PrintClusters() {
  mprintf("CLUSTER: %u clusters.\n", clusters.size());
  for (std::list<clusterNode>::iterator C = clusters.begin(); C != clusters.end(); C++) {
    mprintf("\t%8i : ",(*C).num);
    for (std::list<int>::iterator fnum = (*C).frameList.begin();
                                  fnum != (*C).frameList.end();
                                  fnum++)
    {
      mprintf("%i,",*fnum);
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

// ============================================================================ 
// CONSTRUCTOR
Clustering::Clustering() {
  //fprintf(stderr,"Clustering Con\n");
  useMass=false;
  Linkage = AVERAGELINK;
} 

// DESTRUCTOR
Clustering::~Clustering() { 
}

/* Clustering::init()
 * Expected call: cluster <mask>  
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Clustering::init() {
  char *mask0;

  // Get keywords
  useMass = A->hasKey("mass");
  targetNclusters = A->getKeyInt("clusters",-1);
  epsilon = A->getKeyDouble("epsilon",-1.0);
  if (A->hasKey("linkage")) Linkage=SINGLELINK;
  if (A->hasKey("averagelinkage")) Linkage=SINGLELINK;

  // Get the mask string 
  mask0 = A->getNextMask();
  Mask0.SetMaskString(mask0);

  // Determine finish criteria. If nothing specified default to 10 clusters.
  if (targetNclusters==-1 && epsilon==-1.0)
    targetNclusters = 10;

  mprintf("    CLUSTER: (%s) ",Mask0.maskString);
  if (useMass)
    mprintf(" (mass-weighted)");
  if (targetNclusters != -1)
    mprintf(" looking for %i clusters",targetNclusters);
  if (epsilon != -1.0)
    mprintf(" epsilon is %8.3lf",epsilon);
  mprintf("\n");
  mprintf("            Using hierarchical top-down clustering algorithm,");
  if (Linkage==SINGLELINK)
    mprintf(" single-linkage");
  else if (Linkage==AVERAGELINK)
    mprintf(" average-linkage");
  mprintf(".\n");

  // If epsilon not given make it huge 
  if (epsilon == -1.0) epsilon = DBL_MAX;
  // if target clusters not given make it 1
  if (targetNclusters == -1) targetNclusters=1;

  return 0;
}

/* Clustering::setup()
 * Not important for Clustering, initial pass is only for storing frames.
 */
int Clustering::setup() {
  return 0;  
}

/* Clustering::action()
 * Store current frame as a reference frame.
 */
int Clustering::action() {
  Frame *fCopy;

  fCopy = F->Copy();
  ReferenceFrames.Add(fCopy,P);
  
  return 0;
} 

/* Clustering::calcDistFromRmsd()
 */
int Clustering::calcDistFromRmsd( TriangleMatrix *Distances) {
  // Reference
  AmberParm *RefParm = NULL;
  Frame *RefFrame = NULL;
  Frame *SelectedRef=NULL;
  int lastrefpindex=-1;
  int refatoms = 0;
  // Target
  AmberParm *TgtParm;
  Frame *TgtFrame;
  Frame *SelectedTgt=NULL;
  int lasttgtpindex=-1;
  // Other vars
  double R;
  int current=0;
  int max=0;
  int totalref=0;

  totalref = ReferenceFrames.NumFrames();
  Distances->Setup(totalref);

  max = Distances->Nelements();
  mprintf("  CLUSTER: Calculating RMSDs between each frame (%i total).\n  ",max);

  // Set up progress Bar
  ProgressBar *progress = new ProgressBar(max);

  // LOOP OVER REFERENCE FRAMES
  for (int nref=0; nref < totalref - 1; nref++) {
    progress->Update(current);
    RefParm = ReferenceFrames.GetFrameParm( nref );
    // If the current ref parm not same as last ref parm, reset reference mask
    if (RefParm->pindex != lastrefpindex) {
      if ( Mask0.SetupMask(RefParm,debug) ) {
        mprinterr("Error: Clustering: Could not set up reference mask for %s\n",RefParm->parmName);
        if (SelectedRef!=NULL) delete SelectedRef;
        if (SelectedTgt!=NULL) delete SelectedTgt;
        return 1;
      }
      refatoms = Mask0.Nselected;
      if ( SelectedRef!=NULL ) delete SelectedRef;
      SelectedRef = new Frame(&Mask0, RefParm->mass);
      lastrefpindex = RefParm->pindex;
    }
    // Get the current reference frame
    RefFrame = ReferenceFrames.GetFrame( nref );

    // LOOP OVER TARGET FRAMES
    for (int nframe=nref+1; nframe < totalref; nframe++) {
      TgtParm = ReferenceFrames.GetFrameParm( nframe );
      // If the current frame parm not same as last frame parm, reset frame mask
      if (TgtParm->pindex != lasttgtpindex) {
        if ( Mask0.SetupMask(TgtParm,debug) ) {
          mprinterr("Error: Clustering: Could not set up target mask for %s\n",TgtParm->parmName);
          if (SelectedRef!=NULL) delete SelectedRef;
          if (SelectedTgt!=NULL) delete SelectedTgt;
          return 1;
        }
        // Check that num atoms in mask are the same
        if (Mask0.Nselected != refatoms) {
          mprinterr(
            "Error: Clustering: Num atoms in target mask (%i) != num atoms in ref mask (%i)\n",
                    Mask0.Nselected, refatoms);
          if (SelectedRef!=NULL) delete SelectedRef;
          if (SelectedTgt!=NULL) delete SelectedTgt;
          return 1;
        }
        if ( SelectedTgt!=NULL ) delete SelectedTgt;
        SelectedTgt = new Frame(&Mask0, TgtParm->mass);
        lasttgtpindex = TgtParm->pindex;
      }
      // Get the current target frame
      TgtFrame = ReferenceFrames.GetFrame( nframe );

      // Set selected reference atoms - always done since RMS fit modifies SelectedRef
      SelectedRef->SetFrameCoordsFromMask(RefFrame->X, &Mask0);
      // Set selected target atoms
      SelectedTgt->SetFrameCoordsFromMask(TgtFrame->X, &Mask0);

      // Perform RMS calculation
      R = SelectedTgt->RMSD(SelectedRef, useMass);

      Distances->AddElement( R );
      // DEBUG
      //mprinterr("%12i %12i %12.4lf\n",nref,nframe,R);
      current++;
    } // END loop over target frames
  } // END loop over reference frames
  progress->Update(max);
  delete progress;

  if (SelectedRef!=NULL) delete SelectedRef;
  if (SelectedTgt!=NULL) delete SelectedTgt;

  
  return 0;
}

/* Clustering::ClusterHierAgglo()
 * Cluster using a hierarchical agglomerative (bottom-up) approach. All frames
 * start in their own cluster. The closest two clusters are merged, and 
 * distances between the newly merged cluster and all remaining clusters are
 * recalculated according to one of the following metrics:
 *   single-linkage : The minimum distance between frames in clusters are used.
 *   average-linkage: The average distance between frames in clusters are used.
 */
int Clustering::ClusterHierAgglo( TriangleMatrix *FrameDistances) {
  TriangleMatrix *ClusterDistances;
  ClusterList CList;
  std::list<int> frames;
  bool clusteringComplete = false;
  int C1 = 0;
  int C2 = 0;
  std::list<ClusterList::clusterNode>::iterator C1_it;
  std::list<ClusterList::clusterNode>::iterator C2_it;
  int maxClusters = FrameDistances->Nrows();
  int iterations = 0;
  double min;

  // Build initial clusters.
  for (int cluster = 0; cluster < maxClusters; cluster++) {
    frames.assign(1,cluster);
    CList.AddCluster(&frames,cluster);
  }
  // Build initial cluster distance matrix.
  ClusterDistances = FrameDistances->Copy();

  // DEBUG
  if (debug>1) CList.PrintClusters();

  while (!clusteringComplete) {
    // Find the minimum distance between clusters
    min = ClusterDistances->FindMin(&C1, &C2);
    if (debug>0) mprintf("\tMinimum found between clusters %i and %i (%lf)\n",C1,C2,min);
    // If the minimum distance is greater than epsilon we are done
    if (min > epsilon) {
      mprintf("Minimum distance is greater than epsilon (%lf) clustering complete.\n",epsilon);
      break;
    }

    // Find the clusters in the cluster list
    C1_it = CList.GetCluster(C1);
    C2_it = CList.GetCluster(C2);
    // Sanity check
    //if (C1_it==CList.clusters.end() || C2_it==CList.clusters.end()) return 

    // Merge the closest clusters
    CList.Merge(C1_it,C2_it);
    // DEBUG
    if (debug>1) {
      mprintf("\nAFTER MERGE of %i and %i:\n",C1,C2);
      CList.PrintClusters();
    }
    // Remove all distances having to do with C2
    ClusterDistances->Ignore(C2);

    // If the target number of clusters is reached we are done
    if ((int)CList.clusters.size() <= targetNclusters) {
      mprintf("Target # of clusters (%i) met (%u), clustering complete.\n",targetNclusters,
              CList.clusters.size());
      break;
    } 

    // Recalculate distances between C1 and all other clusters
    if (Linkage==AVERAGELINK)
      CList.calcAvgDist(C1_it, FrameDistances, ClusterDistances);
    else if (Linkage==SINGLELINK)
      CList.calcMinDist(C1_it, FrameDistances, ClusterDistances);
   
    if (debug>2) { 
      mprintf("NEW CLUSTER DISTANCES:\n");
      ClusterDistances->PrintElements();
    }
 
    iterations++;

    // Check if clustering is complete.
    if (CList.clusters.size() == 1) clusteringComplete = true; // Sanity check
  }

  // DEBUG
  mprintf("\nFINAL CLUSTERS after %i iterations:\n",iterations);
  CList.PrintClusters();

  delete ClusterDistances;
  return 0;
}

/* Clustering::print()
 * This is where the clustering is actually performed. First the distances
 * between each frame are calculated. Then the clustering routine is called.
 */
void Clustering::print() {
  TriangleMatrix Distances;

  calcDistFromRmsd( &Distances );

  // DEBUG
  if (debug>1) {
    mprintf("INTIAL FRAME DISTANCES:\n");
    Distances.PrintElements();
  }

  // Cluster
  ClusterHierAgglo( &Distances);
}

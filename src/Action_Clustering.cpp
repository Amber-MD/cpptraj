/* Action: Clustering
 */
#include "Action_Clustering.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#include <cfloat>

// CONSTRUCTOR
Clustering::Clustering() {
  //fprintf(stderr,"Clustering Con\n");
  useMass=false;
  Linkage = AVERAGELINK;
  cnumvtime=NULL;
} 

// DESTRUCTOR
Clustering::~Clustering() { 
}

/* Clustering::init()
 * Expected call: cluster [<mask>] [mass] [clusters <n>] [epsilon <e>] [out <cnumvtime>]
 *                        [ linkage | averagelinkage ]  
 *                        [summary <summaryfile>] 
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Clustering::init() {
  char *mask0,*cnumvtimefile;

  // Get keywords
  useMass = A->hasKey("mass");
  targetNclusters = A->getKeyInt("clusters",-1);
  epsilon = A->getKeyDouble("epsilon",-1.0);
  if (A->hasKey("linkage")) Linkage=SINGLELINK;
  if (A->hasKey("averagelinkage")) Linkage=AVERAGELINK;
  if (A->hasKey("complete")) Linkage=COMPLETELINK;
  cnumvtimefile = A->getKeyString("out",NULL);
  summaryfile = A->getKeyString("summary",NULL);

  // Get the mask string 
  mask0 = A->getNextMask();
  Mask0.SetMaskString(mask0);

  // Dataset to store cluster number v time
  cnumvtime = DSL->Add(INT,A->getNextString(),"Cnum");
  if (cnumvtime==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(cnumvtimefile,cnumvtime);

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
  else if (Linkage==COMPLETELINK)
    mprintf(" complete-linkage");
  mprintf(".\n");
  if (summaryfile!=NULL)
    mprintf("            Summary of cluster results will be written to %s\n",summaryfile);

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
int Clustering::ClusterHierAgglo( TriangleMatrix *FrameDistances, ClusterList *CList) {
  TriangleMatrix *ClusterDistances;
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
    CList->AddCluster(&frames,cluster);
  }
  // Build initial cluster distance matrix.
  ClusterDistances = FrameDistances->Copy();

  // DEBUG
  if (debug>1) CList->PrintClusters();

  while (!clusteringComplete) {
    // Find the minimum distance between clusters
    min = ClusterDistances->FindMin(&C1, &C2);
    if (debug>0) mprintf("\tMinimum found between clusters %i and %i (%lf)\n",C1,C2,min);
    // If the minimum distance is greater than epsilon we are done
    if (min > epsilon) {
      mprintf("\tMinimum distance is greater than epsilon (%lf), clustering complete.\n",
              epsilon);
      break;
    }

    // Find the clusters in the cluster list
    C1_it = CList->GetCluster(C1);
    C2_it = CList->GetCluster(C2);
    // Sanity check
    //if (C1_it==CList.clusters.end() || C2_it==CList.clusters.end()) return 

    // Merge the closest clusters
    CList->Merge(C1_it,C2_it);
    // DEBUG
    if (debug>1) {
      mprintf("\nAFTER MERGE of %i and %i:\n",C1,C2);
      CList->PrintClusters();
    }
    // Remove all distances having to do with C2
    ClusterDistances->Ignore(C2);
  
    // Recalculate distances between C1 and all other clusters
    if (Linkage==AVERAGELINK)
      CList->calcAvgDist(C1_it, FrameDistances, ClusterDistances);
    else if (Linkage==SINGLELINK)
      CList->calcMinDist(C1_it, FrameDistances, ClusterDistances);
    else if (Linkage==COMPLETELINK)
      CList->calcMaxDist(C1_it, FrameDistances, ClusterDistances);
   
    if (debug>2) { 
      mprintf("NEW CLUSTER DISTANCES:\n");
      ClusterDistances->PrintElements();
    }
 
    // Check if clustering is complete.
    // If the target number of clusters is reached we are done
    if ((int)CList->clusters.size() <= targetNclusters) {
      mprintf("\tTarget # of clusters (%i) met (%u), clustering complete.\n",targetNclusters,
              CList->clusters.size());
      break;
    } 
    if (CList->clusters.size() == 1) clusteringComplete = true; // Sanity check
    iterations++;
  }
  mprintf("\tCLUSTER: Completed after %i iterations, %u clusters.\n",iterations,
          CList->clusters.size());


  delete ClusterDistances;
  return 0;
}


/* Clustering::print()
 * This is where the clustering is actually performed. First the distances
 * between each frame are calculated. Then the clustering routine is called.
 */
void Clustering::print() {
  TriangleMatrix Distances;
  ClusterList CList;

  calcDistFromRmsd( &Distances );

  // DEBUG
  if (debug>1) {
    mprintf("INTIAL FRAME DISTANCES:\n");
    Distances.PrintElements();
  }

  // Cluster
  ClusterHierAgglo( &Distances, &CList);

  // Sort clusters and renumber
  CList.Renumber();

  // DEBUG
  if (debug>0) {
    mprintf("\nFINAL CLUSTERS:\n");
    CList.PrintClusters();
  }

  // Print a summary of clusters
  if (summaryfile!=NULL)
    CList.Summary(summaryfile);

  // Create cluster v time data from clusters.
  CList.Cnumvtime( cnumvtime );

}

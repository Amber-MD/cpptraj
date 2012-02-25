/* Action: Clustering
 */
#include "Action_Clustering.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#include "TrajectoryFile.h"
#include <cfloat>
#include <cstdio> // sprintf
#include <cstring> // strlen

// CONSTRUCTOR
Clustering::Clustering() {
  //fprintf(stderr,"Clustering Con\n");
  useMass=false;
  Linkage = ClusterList::AVERAGELINK;
  epsilon=-1.0;
  targetNclusters=-1;
  cnumvtime=NULL;
  summaryfile=NULL;
  halffile=NULL;
  clusterfile=NULL;
  clusterfmt=UNKNOWN_FORMAT;
  singlerepfile=NULL;
  singlerepfmt=UNKNOWN_FORMAT;
  repfile=NULL;
  repfmt=UNKNOWN_FORMAT;
  clusterinfo=NULL;
  nofitrms = false;
  grace_color = false;
  load_pair = true;
  cluster_dataset=NULL;
} 

const char Clustering::PAIRDISTFILE[16]="CpptrajPairDist";

// Clustering::init()
/** Expected call: cluster [<mask>] [mass] [clusters <n>] [epsilon <e>] [out <cnumvtime>]\n 
 *                        [ linkage | averagelinkage | complete ] [gracecolor] [noload] [nofit]\n
 *                        [summary <summaryfile>] [summaryhalf <halffile>] [info <infofile>]\n
 *                        [ clusterout <trajfileprefix> [clusterfmt <trajformat>] ]\n
 *                        [ singlerepout <trajfilename> [singlerepfmt <trajformat>] ]\n
 *                        [ repout <repprefix> [repfmt <repfmt>] ]\n
 */
// Dataset name will be the last arg checked for. Check order is:
//    1) Keywords
//    2) Masks
//    3) Dataset name
int Clustering::init() {
  char *mask0,*cnumvtimefile,*clusterformat,*singlerepformat,*repformat;
  char *dsetname;
  // Get keywords
  useMass = actionArgs.hasKey("mass");
  targetNclusters = actionArgs.getKeyInt("clusters",-1);
  epsilon = actionArgs.getKeyDouble("epsilon",-1.0);
  if (actionArgs.hasKey("linkage")) Linkage=ClusterList::SINGLELINK;
  if (actionArgs.hasKey("averagelinkage")) Linkage=ClusterList::AVERAGELINK;
  if (actionArgs.hasKey("complete")) Linkage=ClusterList::COMPLETELINK;
  cnumvtimefile = actionArgs.getKeyString("out",NULL);
  clusterinfo = actionArgs.getKeyString("info",NULL);
  summaryfile = actionArgs.getKeyString("summary",NULL);
  halffile = actionArgs.getKeyString("summaryhalf",NULL);
  if (actionArgs.hasKey("nofit")) nofitrms=true;
  if (actionArgs.hasKey("gracecolor")) grace_color=true;
  if (actionArgs.hasKey("noload")) load_pair=false;
  if ((dsetname = actionArgs.getKeyString("data",NULL))!=NULL) {
    // Attempt to get dataset from datasetlist
    cluster_dataset = DSL->Get( dsetname );
    if (cluster_dataset == NULL) {
      mprinterr("Error: cluster: dataset %s not found.\n",dsetname);
      return 1;
    }
  }
  // Output trajectory stuff
  clusterfile = actionArgs.getKeyString("clusterout",NULL);
  clusterformat = actionArgs.getKeyString("clusterfmt",NULL);
  singlerepfile = actionArgs.getKeyString("singlerepout",NULL);
  singlerepformat = actionArgs.getKeyString("singlerepfmt",NULL);
  repfile = actionArgs.getKeyString("repout",NULL);
  repformat = actionArgs.getKeyString("repfmt",NULL);
  // Figure out trajectory formats
  if (clusterfile!=NULL) {
    clusterfmt = GetFmtFromArg(clusterformat,AMBERTRAJ);
  }
  if (singlerepfile!=NULL) {
    singlerepfmt = GetFmtFromArg(singlerepformat,AMBERTRAJ);
  }
  if (repfile!=NULL) {
    repfmt = GetFmtFromArg(repformat,AMBERTRAJ);
  }
  // Get the mask string 
  mask0 = actionArgs.getNextMask();
  Mask0.SetMaskString(mask0);

  // Dataset to store cluster number v time
  cnumvtime = DSL->Add(INT,actionArgs.getNextString(),"Cnum");
  if (cnumvtime==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(cnumvtimefile,cnumvtime);

  // Determine finish criteria. If nothing specified default to 10 clusters.
  if (targetNclusters==-1 && epsilon==-1.0)
    targetNclusters = 10;

  mprintf("    CLUSTER:");
  if (dsetname!=NULL) {
    mprintf(" On dataset %s",dsetname);
  } else {
    mprintf(" Using RMSD (mask [%s])",Mask0.MaskString());
    if (useMass)
      mprintf(", mass-weighted");
    if (nofitrms)
      mprintf(", no fitting");
    else
      mprintf(" best fit");
  }
  if (targetNclusters != -1)
    mprintf(" looking for %i clusters",targetNclusters);
  if (epsilon != -1.0)
    mprintf(" epsilon is %8.3lf",epsilon);
  mprintf("\n");
  mprintf("            Using hierarchical bottom-up clustering algorithm,");
  if (Linkage==ClusterList::SINGLELINK)
    mprintf(" single-linkage");
  else if (Linkage==ClusterList::AVERAGELINK)
    mprintf(" average-linkage");
  else if (Linkage==ClusterList::COMPLETELINK)
    mprintf(" complete-linkage");
  mprintf(".\n");
  if (grace_color)
    mprintf("            Grace color instead of cluster number (1-15) will be saved.\n");
  if (load_pair)
    mprintf("            Previously calcd pair distances %s will be used if found.\n",
            PAIRDISTFILE);
  else
    mprintf("            Previously calcd pair distances will be ignored.\n");
  if (clusterinfo!=NULL)
    mprintf("            Cluster information will be written to %s\n",clusterinfo);
  if (summaryfile!=NULL)
    mprintf("            Summary of cluster results will be written to %s\n",summaryfile);
  if (halffile!=NULL)
    mprintf("            Summary comparing first/second half of data for clusters will be written to %s\n",halffile);
  if (clusterfile!=NULL)
    mprintf("            Cluster trajectories will be written to %s, format %s\n",
            clusterfile,File_Format(clusterfmt));
  if (singlerepfile!=NULL)
    mprintf("            Cluster representatives will be written to 1 traj (%s), format %s\n",
            singlerepfile,File_Format(singlerepfmt));
  if (repfile!=NULL) {
    mprintf("            Cluster representatives will be written to separate trajectories,\n");
    mprintf("            prefix (%s), format %s\n",repfile,File_Format(repfmt));
  }
  // If epsilon not given make it huge
  // NOTE: Currently only valid for Hierarchical Agglomerative 
  if (epsilon == -1.0) epsilon = DBL_MAX;
  // if target clusters not given make it 1
  if (targetNclusters == -1) targetNclusters=1;

  return 0;
}

// Clustering::action()
/** Store current frame as a reference frame. Always do this even if
  * not calculating RMSD since we may need to print representative
  * frames etc.
  */
int Clustering::action() {
  Frame *fCopy = currentFrame->FrameCopy();
  ReferenceFrames.AddFrame(fCopy,currentParm);
  
  return 0;
} 

// Clustering::calcDistFromRmsd()
int Clustering::calcDistFromRmsd( TriangleMatrix *Distances) {
  // Reference
  AmberParm *RefParm = NULL;
  Frame *RefFrame = NULL;
  Frame SelectedRef;
  int lastrefpindex=-1;
  int refatoms = 0;
  // Target
  AmberParm *TgtParm;
  Frame *TgtFrame;
  Frame SelectedTgt;
  int lasttgtpindex=-1;
  // Other vars
  double R, U[9], Trans[6];
  int current=0;
  int max=0;
  int totalref=0;

  totalref = ReferenceFrames.NumFrames();
  Distances->Setup(totalref);

  max = Distances->Nelements();
  if (nofitrms)
    mprintf("  CLUSTER: Calculating no-fit RMSDs between each frame (%i total).\n  ",max);
  else
    mprintf("  CLUSTER: Calculating RMSDs with fitting between each frame (%i total).\n",max);

  // Set up progress Bar
  ProgressBar *progress = new ProgressBar(max);

  // LOOP OVER REFERENCE FRAMES
  for (int nref=0; nref < totalref - 1; nref++) {
    progress->Update(current);
    RefParm = ReferenceFrames.GetFrameParm( nref );
    // If the current ref parm not same as last ref parm, reset reference mask
    if (RefParm->pindex != lastrefpindex) {
      if ( RefParm->SetupIntegerMask(Mask0, activeReference)) {
        mprinterr("Error: Clustering: Could not set up reference mask for %s\n",RefParm->parmName);
        return 1;
      }
      refatoms = Mask0.Nselected;
      SelectedRef.SetupFrameFromMask(&Mask0, RefParm->mass);
      lastrefpindex = RefParm->pindex;
    }
    // Get the current reference frame
    RefFrame = ReferenceFrames.GetFrame( nref );

    // LOOP OVER TARGET FRAMES
    for (int nframe=nref+1; nframe < totalref; nframe++) {
      TgtParm = ReferenceFrames.GetFrameParm( nframe );
      // If the current frame parm not same as last frame parm, reset frame mask
      if (TgtParm->pindex != lasttgtpindex) {
        if ( TgtParm->SetupIntegerMask(Mask0, activeReference) ) {
          mprinterr("Error: Clustering: Could not set up target mask for %s\n",TgtParm->parmName);
          return 1;
        }
        // Check that num atoms in mask are the same
        if (Mask0.Nselected != refatoms) {
          mprinterr(
            "Error: Clustering: Num atoms in target mask (%i) != num atoms in ref mask (%i)\n",
                    Mask0.Nselected, refatoms);
          return 1;
        }
        SelectedTgt.SetupFrameFromMask(&Mask0, TgtParm->mass);
        lasttgtpindex = TgtParm->pindex;
      }
      // Get the current target frame
      TgtFrame = ReferenceFrames.GetFrame( nframe );

      // Set selected reference atoms - always done since RMS fit modifies SelectedRef
      SelectedRef.SetFrameCoordsFromMask(RefFrame->X, &Mask0);
      // Set selected target atoms
      SelectedTgt.SetFrameCoordsFromMask(TgtFrame->X, &Mask0);

      // Perform RMS calculation
      if (nofitrms)
        R = SelectedTgt.RMSD(&SelectedRef, useMass);
      else 
        R = SelectedTgt.RMSD(&SelectedRef, U, Trans, useMass);

      Distances->AddElement( R );
      // DEBUG
      //mprinterr("%12i %12i %12.4lf\n",nref,nframe,R);
      current++;
    } // END loop over target frames
  } // END loop over reference frames
  progress->Update(max);
  delete progress;

  return 0;
}

// Clustering::ClusterHierAgglo()
/** Cluster using a hierarchical agglomerative (bottom-up) approach. All frames
  * start in their own cluster. The closest two clusters are merged, and 
  * distances between the newly merged cluster and all remaining clusters are
  * recalculated according to one of the following metrics:
  * - single-linkage  : The minimum distance between frames in clusters are used.
  * - average-linkage : The average distance between frames in clusters are used.
  * - complete-linkage: The maximum distance between frames in clusters are used.
  */
int Clustering::ClusterHierAgglo( TriangleMatrix *FrameDistances, ClusterList *CList) {
  std::list<int> frames;
  bool clusteringComplete = false;
  int iterations = 0;

  mprintf("\tStarting Hierarchical Agglomerative Clustering:\n");
  ProgressBar cluster_progress(-1);

  // Build initial clusters.
  for (int cluster = 0; cluster < FrameDistances->Nrows(); cluster++) {
    frames.assign(1,cluster);
    CList->AddCluster(&frames,cluster);
  }
  // Build initial cluster distance matrix.
  CList->Initialize( FrameDistances );

  // DEBUG
  if (debug>1) CList->PrintClusters();

  // Initial check to see if epsilon is satisfied.
  //clusteringComplete = CList->CheckEpsilon(epsilon);

  while (!clusteringComplete) {
    // Merge 2 closest clusters
    if (CList->MergeClosest(epsilon)) break; 

    // Check if clustering is complete.
    // If the target number of clusters is reached we are done
    if (CList->Nclusters() <= targetNclusters) {
      mprintf("\tTarget # of clusters (%i) met (%u), clustering complete.\n",targetNclusters,
              CList->Nclusters());
      break;
    } 
    if (CList->Nclusters() == 1) clusteringComplete = true; // Sanity check
    cluster_progress.Update( iterations );
    iterations++;
  }
  mprintf("\tCLUSTER: Completed after %i iterations, %u clusters.\n",iterations,
          CList->Nclusters());

  return 0;
}

// Clustering::CreateCnumvtime()
/** Put cluster number vs frame into dataset.  */
void Clustering::CreateCnumvtime( ClusterList *CList ) {
  std::list<int>::iterator E;
  int cnum;

  CList->Begin();
  while (!CList->End()) {
    //mprinterr("Cluster %i:\n",CList->CurrentNum());
    cnum = CList->CurrentNum();
    // If grace colors, return integer in range from 1 to 15 (1 most populated)
    if (grace_color) {
      cnum = cnum + 1;
      if (cnum > 15) cnum = 15;
    } 
    // Loop over all frames in the cluster
    E = CList->CurrentFrameEnd();
    for (std::list<int>::iterator frame = CList->CurrentFrameBegin();
                                  frame != E;
                                  frame++)
    {
      //mprinterr("%i,",*frame);
      cnumvtime->Add( *frame, &cnum );
    }
    //mprinterr("\n");
    CList->NextCluster();
    //break;
  }
}

// Clustering::WriteClusterTraj()
/** Write frames in each cluster to a trajectory file.  */
void Clustering::WriteClusterTraj( ClusterList *CList ) {
  std::list<int>::iterator E;
  std::list<int>::iterator B;
  char *cfilename;
  int cnum,framenum;
  TrajectoryFile *clusterout = NULL;
  AmberParm *clusterparm;
  Frame *clusterframe;

  // Figure out max size of cluster filename
  cnum = CList->Nclusters();
  cnum = (cnum / 10) + 3;
  cfilename = new char[ strlen(clusterfile)+cnum+1 ];

  CList->Begin();
  while (!CList->End()) {
    // Create filename based on cluster number.
    cnum = CList->CurrentNum();
    sprintf(cfilename,"%s.c%i",clusterfile,cnum);
    // Set up trajectory file - use parm from first frame of cluster (pot. dangerous)
    if (clusterout!=NULL) delete clusterout;
    clusterout = new TrajectoryFile;
    B = CList->CurrentFrameBegin();
    clusterparm = ReferenceFrames.GetFrameParm( *B );
    if (clusterout->SetupWrite(cfilename,NULL,clusterparm,clusterfmt)) {
      mprinterr("Error: Clustering::WriteClusterTraj: Could not set up %s for write.\n",
                cfilename);
      delete clusterout;
      delete[] cfilename;
      return;
    }
    //mprinterr("Cluster %i:\n",CList->CurrentNum());
    
    E = CList->CurrentFrameEnd();
    framenum = 0;
    for (std::list<int>::iterator frame = B; frame != E; frame++)
    {
      //mprinterr("%i,",*frame);
      clusterframe = ReferenceFrames.GetFrame( *frame );
      clusterout->WriteFrame(framenum++, clusterparm, *clusterframe);
    }
    // Close traj
    clusterout->EndTraj();
    //mprinterr("\n");
    CList->NextCluster();
    //break;
  }
  delete[] cfilename;
  if (clusterout!=NULL) delete clusterout;
}

// Clustering::WriteSingleRepTraj()
/** Write representative frame of each cluster to a trajectory file.  */
void Clustering::WriteSingleRepTraj( ClusterList *CList ) {
  int framenum, framecounter;
  TrajectoryFile clusterout;
  AmberParm *clusterparm;
  Frame *clusterframe;

  // Find centroid of first cluster in order to set up parm
  // NOTE: This is redundant if the Summary routine has already been called.
  CList->Begin();
  framenum = CList->CurrentCentroid();

  // Set up trajectory file. Use parm from first frame of cluster (pot. dangerous)
  clusterparm = ReferenceFrames.GetFrameParm( framenum );
  if (clusterout.SetupWrite(singlerepfile,NULL,clusterparm,singlerepfmt)) {
    mprinterr("Error: Clustering::WriteSingleRepTraj: Could not set up %s for write.\n",
                singlerepfile);
     return;
  }
  // Write first cluster rep frame
  framecounter=0;
  clusterframe = ReferenceFrames.GetFrame( framenum );
  clusterout.WriteFrame(framecounter++, clusterparm, *clusterframe);

  CList->NextCluster();
  while (!CList->End()) {
    //mprinterr("Cluster %i: ",CList->CurrentNum());
   framenum = CList->CurrentCentroid();
   //mprinterr("%i\n",framenum);
   clusterframe = ReferenceFrames.GetFrame( framenum );
   clusterout.WriteFrame(framecounter++, clusterparm, *clusterframe);
    //mprinterr("\n");
    CList->NextCluster();
    //break;
  }
  // Close traj
  clusterout.EndTraj();
}

// Clustering::WriteRepTraj()
/** Write representative frame of each cluster to a separate trajectory file,
  * repfile.REPNUM.FMT
  */
void Clustering::WriteRepTraj( ClusterList *CList ) {
  int framenum, cnum;
  TrajectoryFile *clusterout = NULL;
  AmberParm *clusterparm;
  Frame *clusterframe;
  char ext[8];
  char *cfilename = NULL;

  // Set output extension for this file format
  SetExtFromFmt(ext,repfmt);

  CList->Begin();
  while (!CList->End()) {
    // Create trajectory file object
    if (clusterout!=NULL) delete clusterout;
    clusterout = new TrajectoryFile();

    // Find centroid of first cluster in order to set up parm
    framenum = CList->CurrentCentroid();

    // Set up trajectory filename for this rep frame
    cnum = (framenum / 10) + 8;
    if (cfilename!=NULL) delete[] cfilename;
    cfilename = new char[ strlen(repfile)+cnum+1 ];
    sprintf(cfilename, "%s.%i%s", repfile, framenum+1, ext);

    // Set up trajectory file. Use parm from first frame of cluster (pot. dangerous)
    clusterparm = ReferenceFrames.GetFrameParm( framenum );
    if (clusterout->SetupWrite(cfilename,NULL,clusterparm,repfmt)) {
      mprinterr("Error: Clustering::WriteRepTraj: Could not set up %s for write.\n",
                cfilename);
       return;
    }

    // Write cluster rep frame
    clusterframe = ReferenceFrames.GetFrame( framenum );
    clusterout->WriteFrame(framenum, clusterparm, *clusterframe);
    // Close traj
    clusterout->EndTraj();

    CList->NextCluster();
  }
}

// Clustering::calcDistFromDataset()
void Clustering::calcDistFromDataset( TriangleMatrix &Distances ) {
  int N = cluster_dataset->Xmax();
  // Since Xmax returns the last point added, add 1
  ++N;
  //mprintf("DEBUG: xmax is %i\n",N);
  Distances.Setup(N);

  int max = Distances.Nelements();
  mprintf("  CLUSTER: Calculating distances using dataset %s (%i total).\n",
          cluster_dataset->Name(),max);

  ProgressBar *progress = new ProgressBar(max);

  // LOOP 
  int current = 0;
  for (int i = 0; i < N-1; i++) {
    progress->Update(current);
    double iVal = cluster_dataset->Dval(i);
    for (int j = i + 1; j < N; j++) {
      double jVal = cluster_dataset->Dval(j);
      // Calculate abs( delta )
      double delta = iVal - jVal;
      if (delta < 0) delta = -delta;
      Distances.AddElement( delta );
      current++;
    }
  }
  progress->Update(max);
  delete progress;

}

// Clustering::print()
/** This is where the clustering is actually performed. First the distances
  * between each frame are calculated. Then the clustering routine is called.
  */
void Clustering::print() {
  TriangleMatrix Distances;
  ClusterList CList;

  // If PAIRDISTFILE exists load pair distances from there
  if (load_pair && fileExists((char*)PAIRDISTFILE)) {
    mprintf("CLUSTER: %s found, loading pairwise distances.\n",PAIRDISTFILE);
    if (Distances.LoadFile((char*)PAIRDISTFILE,ReferenceFrames.NumFrames())) return;
  } else if (cluster_dataset==NULL) {
    // Get RMSDs between frames
    calcDistFromRmsd( &Distances );
    // Save distances
    // NOTE: Only if load_pair?
    Distances.SaveFile((char*)PAIRDISTFILE);
  } else {
    // Get distances from dataset
    calcDistFromDataset( Distances );
    // NOTE: Need to update save to indicate distance type
    Distances.SaveFile((char*)PAIRDISTFILE);
  }

  // DEBUG
  if (debug>1) {
    mprintf("INTIAL FRAME DISTANCES:\n");
    Distances.PrintElements();
  }

  // Cluster
  CList.SetDebug(debug);
  CList.SetLinkage(Linkage);
  ClusterHierAgglo( &Distances, &CList);

  // Sort clusters and renumber; also finds centroids for printing
  // representative frames.
  CList.Renumber();

  // DEBUG
  if (debug>0) {
    mprintf("\nFINAL CLUSTERS:\n");
    CList.PrintClusters();
    mprintf("\nREPRESENTATIVE FRAMES:\n");
    CList.PrintRepFrames();
  }

  // Print ptraj-like cluster info
  if (clusterinfo!=NULL)
    CList.PrintClustersToFile(clusterinfo);

  // Print a summary of clusters
  if (summaryfile!=NULL)
    CList.Summary(summaryfile);

  // Print a summary comparing first half to second half of data for clusters
  if (halffile!=NULL)
    CList.Summary_Half(halffile);

  // Create cluster v time data from clusters.
  CreateCnumvtime( &CList );

  // Write clusters to trajectories
  if (clusterfile!=NULL)
    WriteClusterTraj( &CList ); 

  // Write all representative frames to a single traj
  if (singlerepfile!=NULL)
    WriteSingleRepTraj( &CList );

  // Write all representative frames to separate trajs
  if (repfile!=NULL)
    WriteRepTraj( &CList );
}

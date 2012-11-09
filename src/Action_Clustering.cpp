// Action: Action_Clustering
#include <cfloat> // DBL_MAX
#include "Action_Clustering.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // fileExists, integerToString
#include "ProgressBar.h"
#include "DataSet_integer.h" // For converting cnumvtime
#include "Trajout.h"

// CONSTRUCTOR
Action_Clustering::Action_Clustering() :
  epsilon_(-1.0),
  targetNclusters_(-1),
  sieve_(1),
  cnumvtime_(NULL),
  nofitrms_(false),
  useMass_(false),
  grace_color_(false),
  load_pair_(true),
  cluster_dataset_(NULL),
  Linkage_(ClusterList::AVERAGELINK),
  clusterfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  singlerepfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  reptrajfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  CurrentParm_(0)
{ } 

void Action_Clustering::Help() {
  mprintf("cluster [<mask>] [mass] [clusters <n>] [epsilon <e>] [out <cnumvtime>]\n");
  mprintf("        [ linkage | averagelinkage | complete ] [gracecolor] [noload] [nofit]\n");
  mprintf("        [summary <summaryfile>] [summaryhalf <halffile>] [info <infofile>]\n");
  mprintf("        [ clusterout <trajfileprefix> [clusterfmt <trajformat>] ]\n");
  mprintf("        [ singlerepout <trajfilename> [singlerepfmt <trajformat>] ]\n");
  mprintf("        [ repout <repprefix> [repfmt <repfmt>] ]\n");
  mprintf("        [data <setname>]\n");
  mprintf("\tCluster structures based on RMSD or a given DataFile.\n");
}

const char Action_Clustering::PAIRDISTFILE[16]="CpptrajPairDist";

// Action_Clustering::init()
Action::RetType Action_Clustering::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Get keywords
  useMass_ = actionArgs.hasKey("mass");
  targetNclusters_ = actionArgs.getKeyInt("clusters",-1);
  sieve_ = actionArgs.getKeyInt("sieve",1);
  epsilon_ = actionArgs.getKeyDouble("epsilon",-1.0);
  if (actionArgs.hasKey("linkage")) Linkage_=ClusterList::SINGLELINK;
  if (actionArgs.hasKey("averagelinkage")) Linkage_=ClusterList::AVERAGELINK;
  if (actionArgs.hasKey("complete")) Linkage_=ClusterList::COMPLETELINK;
  ArgList::ConstArg cnumvtimefile = actionArgs.getKeyString("out");
  clusterinfo_ = actionArgs.GetStringKey("info");
  summaryfile_ = actionArgs.GetStringKey("summary");
  halffile_ = actionArgs.GetStringKey("summaryhalf");
  if (actionArgs.hasKey("nofit")) nofitrms_=true;
  if (actionArgs.hasKey("gracecolor")) grace_color_=true;
  if (actionArgs.hasKey("noload")) load_pair_=false;
  std::string dsetname = actionArgs.GetStringKey("data");
  if (!dsetname.empty()) {
    // Attempt to get dataset from datasetlist
    cluster_dataset_ = DSL->GetDataSet( dsetname );
    if (cluster_dataset_ == 0) {
      mprinterr("Error: cluster: dataset %s not found.\n",dsetname.c_str());
      return Action::ERR;
    }
  }
  // Output trajectory stuff
  clusterfile_ = actionArgs.GetStringKey("clusterout");
  clusterfmt_ = TrajectoryFile::GetFormatFromString( actionArgs.GetStringKey("clusterfmt") ); 
  singlerepfile_ = actionArgs.GetStringKey("singlerepout");
  singlerepfmt_ = TrajectoryFile::GetFormatFromString( actionArgs.GetStringKey("singlerepfmt") );
  reptrajfile_ = actionArgs.GetStringKey("repout");
  reptrajfmt_ = TrajectoryFile::GetFormatFromString( actionArgs.GetStringKey("repfmt") );
  // Get the mask string 
  Mask0_.SetMaskString( actionArgs.getNextMask() );

  // Dataset to store cluster number v time
  cnumvtime_ = DSL->Add(DataSet::INT, actionArgs.getNextString(), "Cnum");
  if (cnumvtime_==NULL) return Action::ERR;
  // Add dataset to data file list
  DFL->Add(cnumvtimefile,cnumvtime_);

  // Determine finish criteria. If nothing specified default to 10 clusters.
  if (targetNclusters_==-1 && epsilon_==-1.0)
    targetNclusters_ = 10;

  mprintf("    CLUSTER:");
  if (!dsetname.empty()) {
    mprintf(" On dataset %s",dsetname.c_str());
  } else {
    mprintf(" Using RMSD (mask [%s])",Mask0_.MaskString());
    if (useMass_)
      mprintf(", mass-weighted");
    if (nofitrms_)
      mprintf(", no fitting");
    else
      mprintf(" best fit");
  }
  if (targetNclusters_ != -1)
    mprintf(" looking for %i clusters",targetNclusters_);
  if (epsilon_ != -1.0)
    mprintf(" epsilon is %.3lf",epsilon_);
  mprintf("\n");
  mprintf("\tUsing hierarchical bottom-up clustering algorithm,");
  if (Linkage_==ClusterList::SINGLELINK)
    mprintf(" single-linkage");
  else if (Linkage_==ClusterList::AVERAGELINK)
    mprintf(" average-linkage");
  else if (Linkage_==ClusterList::COMPLETELINK)
    mprintf(" complete-linkage");
  mprintf(".\n");
  if (grace_color_)
    mprintf("\tGrace color instead of cluster number (1-15) will be saved.\n");
  if (load_pair_)
    mprintf("\tPreviously calcd pair distances %s will be used if found.\n",
            PAIRDISTFILE);
  else
    mprintf("\tPreviously calcd pair distances will be ignored.\n");
  if (!clusterinfo_.empty())
    mprintf("\tCluster information will be written to %s\n",clusterinfo_.c_str());
  if (!summaryfile_.empty())
    mprintf("\tSummary of cluster results will be written to %s\n",summaryfile_.c_str());
  if (!halffile_.empty())
    mprintf("\tSummary comparing first/second half of data for clusters will be written to %s\n",
            halffile_.c_str());
  if (!clusterfile_.empty())
    mprintf("\tCluster trajectories will be written to %s, format %s\n",
            clusterfile_.c_str(), TrajectoryFile::FormatString(clusterfmt_));
  if (!singlerepfile_.empty())
    mprintf("\tCluster representatives will be written to 1 traj (%s), format %s\n",
            singlerepfile_.c_str(), TrajectoryFile::FormatString(singlerepfmt_));
  if (!reptrajfile_.empty()) {
    mprintf("\tCluster representatives will be written to separate trajectories,\n");
    mprintf("\t\tprefix (%s), format %s\n",reptrajfile_.c_str(), 
            TrajectoryFile::FormatString(reptrajfmt_));
  }
  // If epsilon not given make it huge
  // NOTE: Currently only valid for Hierarchical Agglomerative 
  if (epsilon_ == -1.0) epsilon_ = DBL_MAX;
  // if target clusters not given make it 1
  if (targetNclusters_ == -1) targetNclusters_=1;

  return Action::OK;
}

Action::RetType Action_Clustering::Setup(Topology* currentParm, Topology** parmAddress) {
  CurrentParm_ = currentParm;
  return Action::OK;
}

// Action_Clustering::action()
/** Store current frame as a reference frame. Always do this even if
  * not calculating RMSD since we may need to print representative
  * frames etc.
  */
Action::RetType Action_Clustering::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  Frame *fCopy = currentFrame->FrameCopy();
  ReferenceFrames_.AddFrame(fCopy,CurrentParm_);
  
  return Action::OK;
} 

// Action_Clustering::print()
/** This is where the clustering is actually performed. First the distances
  * between each frame are calculated. Then the clustering routine is called.
  */
// TODO: Need to update save to indicate distance type
// NOTE: Should distances be saved only if load_pair?
void Action_Clustering::Print() {
  TriangleMatrix Distances;
  ClusterList CList;

  mprintf("    CLUSTER:");
  // Default: 0 - Get RMSDs from frames.
  //          1 - If PAIRDISTFILE exists load pair distances from there.
  //          2 - If DataSet was specified get pair distances from that.
  // Calculated distances will be saved if not loaded from file.
  int pairdist_mode = 0; 
  if (load_pair_ && fileExists(PAIRDISTFILE))
    pairdist_mode = 1;
  if (cluster_dataset_ != NULL)
    pairdist_mode = 2;

  if (pairdist_mode == 2) {  // Get distances from dataset.
    calcDistFromDataset( Distances );
    Distances.SaveFile( PAIRDISTFILE );
  }
  if (pairdist_mode == 1) {  // Get distances from file.
    mprintf(" %s found, loading pairwise distances.\n",PAIRDISTFILE);
    if (Distances.LoadFile(PAIRDISTFILE,ReferenceFrames_.NumFrames())) {
      mprintf("\tLoading pairwise distances failed - regenerating from frames.\n");
      pairdist_mode = 0;
    }
  } 
  if (pairdist_mode == 0) { // Get RMSDs between frames
    calcDistFromRmsd( Distances );
    Distances.SaveFile( PAIRDISTFILE );
  } 

  // DEBUG
  if (debug_>1) {
    mprintf("INTIAL FRAME DISTANCES:\n");
    Distances.PrintElements();
  }

  // Cluster
  CList.SetDebug(debug_);
  CList.SetLinkage(Linkage_);
  ClusterHierAgglo( Distances, CList);

  // Sort clusters and renumber; also finds centroids for printing
  // representative frames.
  CList.Renumber();

  // DEBUG
  if (debug_>0) {
    mprintf("\nFINAL CLUSTERS:\n");
    CList.PrintClusters();
    mprintf("\nREPRESENTATIVE FRAMES:\n");
    CList.PrintRepFrames();
  }

  // Print ptraj-like cluster info
  if (!clusterinfo_.empty())
    CList.PrintClustersToFile(clusterinfo_);

  // Print a summary of clusters
  if (!summaryfile_.empty())
    CList.Summary(summaryfile_);

  // Print a summary comparing first half to second half of data for clusters
  if (!halffile_.empty())
    CList.Summary_Half(halffile_);

  // Create cluster v time data from clusters.
  CreateCnumvtime( CList );

  // Write clusters to trajectories
  if (!clusterfile_.empty())
    WriteClusterTraj( CList ); 

  // Write all representative frames to a single traj
  if (!singlerepfile_.empty())
    WriteSingleRepTraj( CList );

  // Write all representative frames to separate trajs
  if (!reptrajfile_.empty())
    WriteRepTraj( CList );
}

// -----------------------------------------------------------------------------
// Action_Clustering::calcDistFromRmsd()
int Action_Clustering::calcDistFromRmsd( TriangleMatrix& Distances) {
  // Reference
  Topology *RefParm = NULL;
  Frame *RefFrame = NULL;
  Frame SelectedRef;
  int lastrefpindex=-1;
  int refatoms = 0;
  // Target
  Topology *TgtParm;
  Frame *TgtFrame;
  Frame SelectedTgt;
  int lasttgtpindex=-1;
  // Other vars
  double R, U[9], Trans[6];
  int current=0;
  int max=0;
  int totalref=0;

  totalref = ReferenceFrames_.NumFrames();
  Distances.Setup(totalref);

  max = Distances.Nelements();
  if (nofitrms_)
    mprintf(" Calculating no-fit RMSDs between each frame (%i total).\n  ",max);
  else
    mprintf(" Calculating RMSDs with fitting between each frame (%i total).\n",max);

  // Set up progress Bar
  ProgressBar progress(max);

  // LOOP OVER REFERENCE FRAMES
  for (int nref=0; nref < totalref - 1; nref++) {
    progress.Update(current);
    RefParm = ReferenceFrames_.GetFrameParm( nref );
    // If the current ref parm not same as last ref parm, reset reference mask
    if (RefParm->Pindex() != lastrefpindex) {
      if ( RefParm->SetupIntegerMask(Mask0_)) {
        mprinterr("Error: Clustering: Could not set up reference mask for %s\n",RefParm->c_str());
        return 1;
      }
      refatoms = Mask0_.Nselected();
      // NOTE: This copies in correct masses according to Mask0
      SelectedRef.SetupFrameFromMask(Mask0_, RefParm->Atoms());
      lastrefpindex = RefParm->Pindex();
    }
    // Get the current reference frame
    RefFrame = ReferenceFrames_.GetFrame( nref );
    // Set the selected atoms from the reference frame
    SelectedRef.SetCoordinates(*RefFrame, Mask0_);
    // If fitting, pre-center reference frame
    if (!nofitrms_)
      SelectedRef.CenterReference(Trans+3, useMass_);

    // LOOP OVER TARGET FRAMES
    for (int nframe=nref+1; nframe < totalref; nframe++) {
      TgtParm = ReferenceFrames_.GetFrameParm( nframe );
      // If the current frame parm not same as last frame parm, reset frame mask
      if (TgtParm->Pindex() != lasttgtpindex) {
        if ( TgtParm->SetupIntegerMask(Mask0_) ) {
          mprinterr("Error: Clustering: Could not set up target mask for %s\n",TgtParm->c_str());
          return 1;
        }
        // Check that num atoms in mask are the same
        if (Mask0_.Nselected() != refatoms) {
          mprintf("Warning: Clustering RMS: Num atoms in frame %i (%i) != num atoms \n",
                    nframe+1, Mask0_.Nselected());
          mprintf("Warning:   in frame %i (%i). Assigning an RMS of 10000.\n",
                   nref+1, refatoms);
          R = 10000;
          Distances.AddElement( R );
          continue;
        }
        // NOTE: This copies in correct masses according to Mask0
        SelectedTgt.SetupFrameFromMask(Mask0_, TgtParm->Atoms());
        lasttgtpindex = TgtParm->Pindex();
      }
      // Get the current target frame
      TgtFrame = ReferenceFrames_.GetFrame( nframe );

      // Set the selected atoms from the target frame
      SelectedTgt.SetCoordinates(*TgtFrame, Mask0_);

      // Perform RMS calculation
      if (nofitrms_)
        R = SelectedTgt.RMSD(SelectedRef, useMass_);
      else 
        R = SelectedTgt.RMSD_CenteredRef(SelectedRef, U, Trans, useMass_);

      Distances.AddElement( R );
      // DEBUG
      //mprinterr("%12i %12i %12.4lf\n",nref,nframe,R);
      ++current;
    } // END loop over target frames
  } // END loop over reference frames

  return 0;
}

// Action_Clustering::ClusterHierAgglo()
/** Cluster using a hierarchical agglomerative (bottom-up) approach. All frames
  * start in their own cluster. The closest two clusters are merged, and 
  * distances between the newly merged cluster and all remaining clusters are
  * recalculated according to one of the following metrics:
  * - single-linkage  : The minimum distance between frames in clusters are used.
  * - average-linkage : The average distance between frames in clusters are used.
  * - complete-linkage: The maximum distance between frames in clusters are used.
  */
int Action_Clustering::ClusterHierAgglo( TriangleMatrix& FrameDistances, 
                                  ClusterList& CList) 
{
  std::list<int> frames;
  bool clusteringComplete = false;
  int iterations = 0;

  mprintf("\tStarting Hierarchical Agglomerative Clustering:\n");
  ProgressBar cluster_progress(-1);

  // Build initial clusters.
  for (int cluster = 0; cluster < FrameDistances.Nrows(); cluster++) {
    frames.assign(1,cluster);
    CList.AddCluster(frames, cluster);
  }
  mprintf("\t%i initial clusters.\n", CList.Nclusters());
  // Build initial cluster distance matrix.
  CList.Initialize( &FrameDistances );

  // DEBUG
  if (debug_>1) CList.PrintClusters();

  // Initial check to see if epsilon is satisfied.
  //clusteringComplete = CList->CheckEpsilon(epsilon_);

  while (!clusteringComplete) {
    // Merge 2 closest clusters
    if (CList.MergeClosest(epsilon_)) break; 

    // Check if clustering is complete.
    // If the target number of clusters is reached we are done
    if (CList.Nclusters() <= targetNclusters_) {
      mprintf("\n\tTarget # of clusters (%i) met (%u), clustering complete.\n",targetNclusters_,
              CList.Nclusters());
      break;
    } 
    if (CList.Nclusters() == 1) clusteringComplete = true; // Sanity check
    cluster_progress.Update( iterations );
    ++iterations;
  }
  mprintf("\tCompleted after %i iterations, %u clusters.\n",iterations,
          CList.Nclusters());

  return 0;
}

// Action_Clustering::CreateCnumvtime()
/** Put cluster number vs frame into dataset.  */
void Action_Clustering::CreateCnumvtime( ClusterList &CList ) {
  // FIXME:
  // Cast generic DataSet for cnumvtime back to integer dataset to 
  // access specific integer dataset functions for resizing and []
  // operator. Should this eventually be generic to all atomic DataSets? 
  DataSet_integer* cnum_temp = (DataSet_integer*)cnumvtime_;
  cnum_temp->Resize( CList.MaxFrames() );

  for (ClusterList::cluster_iterator C = CList.begincluster();
                                     C != CList.endcluster(); C++)
  {
    //mprinterr("Cluster %i:\n",CList->CurrentNum());
    int cnum = (*C).Num();
    // If grace colors, return integer in range from 1 to 15 (1 most populated)
    if (grace_color_) {
      cnum = cnum + 1;
      if (cnum > 15) cnum = 15;
    } 
    // Loop over all frames in the cluster
    for (ClusterNode::frame_iterator frame = (*C).beginframe();
                                     frame != (*C).endframe(); frame++)
    {
      //mprinterr("%i,",*frame);
      (*cnum_temp)[ *frame ] = cnum;
    }
    //mprinterr("\n");
    //break;
  }
}

// Action_Clustering::WriteClusterTraj()
/** Write frames in each cluster to a trajectory file.  */
void Action_Clustering::WriteClusterTraj( ClusterList &CList ) {
  // Loop over all clusters
  for (ClusterList::cluster_iterator C = CList.begincluster();
                                     C != CList.endcluster(); C++)
  {
    // Create filename based on cluster number.
    int cnum = (*C).Num();
    std::string cfilename =  clusterfile_ + ".c" + integerToString( cnum );
    // Set up trajectory file 
    // Use parm from first frame of cluster (pot. dangerous)
    Trajout *clusterout = new Trajout;
    ClusterNode::frame_iterator frame = (*C).beginframe();
    Topology *clusterparm = ReferenceFrames_.GetFrameParm( *frame );
    if (clusterout->SetupTrajWrite(cfilename, NULL, clusterparm, clusterfmt_)) 
    {
      mprinterr("Error: Clustering::WriteClusterTraj: Could not set up %s for write.\n",
                cfilename.c_str());
      delete clusterout;
      return;
    }
    //mprinterr("Cluster %i:\n",CList->CurrentNum());
    // Loop over all frames in cluster
    int framenum = 0;
    for (; frame != (*C).endframe(); frame++) {
      //mprinterr("%i,",*frame);
      Frame *clusterframe = ReferenceFrames_.GetFrame( *frame );
      clusterout->WriteFrame(framenum++, clusterparm, *clusterframe);
    }
    // Close traj
    clusterout->EndTraj();
    //mprinterr("\n");
    //break;
    delete clusterout;
  }
}

// Action_Clustering::WriteSingleRepTraj()
/** Write representative frame of each cluster to a trajectory file.  */
void Action_Clustering::WriteSingleRepTraj( ClusterList &CList ) {
  Trajout clusterout;

  // Find centroid of first cluster in order to set up parm
  // NOTE: This is redundant if the Summary routine has already been called.
  ClusterList::cluster_iterator cluster = CList.begincluster();
  int framenum = (*cluster).Centroid();

  // Set up trajectory file. Use parm from first frame of cluster (pot. dangerous)
  Topology *clusterparm = ReferenceFrames_.GetFrameParm( framenum );
  if (clusterout.SetupTrajWrite(singlerepfile_, NULL, clusterparm, singlerepfmt_)) 
  {
    mprinterr("Error: Clustering::WriteSingleRepTraj: Could not set up %s for write.\n",
                singlerepfile_.c_str());
     return;
  }
  // Write first cluster rep frame
  int framecounter = 0;
  Frame *clusterframe = ReferenceFrames_.GetFrame( framenum );
  clusterout.WriteFrame(framecounter++, clusterparm, *clusterframe);

  ++cluster;
  for (; cluster != CList.endcluster(); cluster++) {
    //mprinterr("Cluster %i: ",CList->CurrentNum());
   framenum = (*cluster).Centroid();
   //mprinterr("%i\n",framenum);
   clusterframe = ReferenceFrames_.GetFrame( framenum );
   clusterout.WriteFrame(framecounter++, clusterparm, *clusterframe);
    //mprinterr("\n");
    //break;
  }
  // Close traj
  clusterout.EndTraj();
}

// Action_Clustering::WriteRepTraj()
/** Write representative frame of each cluster to a separate trajectory file,
  * repfile.REPNUM.FMT
  */
void Action_Clustering::WriteRepTraj( ClusterList &CList ) {
  // Create trajectory file object
  std::string tmpExt = TrajectoryFile::GetExtensionForType(reptrajfmt_);

  for (ClusterList::cluster_iterator C = CList.begincluster();
                                     C != CList.endcluster(); ++C)
  {
    Trajout *clusterout = new Trajout();
    // Find centroid of first cluster in order to set up parm
    int framenum = (*C).Centroid();
    // Create filename
    std::string cfilename = reptrajfile_ + "." + integerToString(framenum+1) + tmpExt;
    // Set up trajectory file. 
    // Use parm from first frame of cluster (pot. dangerous)
    Topology *clusterparm = ReferenceFrames_.GetFrameParm( framenum );
    if (clusterout->SetupTrajWrite(cfilename, NULL, clusterparm, reptrajfmt_)) 
    {
      mprinterr("Error: Clustering::WriteRepTraj: Could not set up %s for write.\n",
                reptrajfile_.c_str());
       delete clusterout;
       return;
    }
    // Write cluster rep frame
    Frame *clusterframe = ReferenceFrames_.GetFrame( framenum );
    clusterout->WriteFrame(framenum, clusterparm, *clusterframe);
    // Close traj
    clusterout->EndTraj();
    delete clusterout;
  }
}

// Action_Clustering::calcDistFromDataset()
void Action_Clustering::calcDistFromDataset( TriangleMatrix &Distances ) {
  int N = cluster_dataset_->Size();
  //mprintf("DEBUG: xmax is %i\n",N);
  Distances.Setup(N);

  int max = Distances.Nelements();
  mprintf(" Calculating distances using dataset %s (%i total).\n",
          cluster_dataset_->Legend().c_str(),max);

  ProgressBar progress(max);
  // LOOP 
  int current = 0;
  for (int i = 0; i < N-1; i++) {
    progress.Update(current);
    double iVal = cluster_dataset_->Dval(i);
    for (int j = i + 1; j < N; j++) {
      double jVal = cluster_dataset_->Dval(j);
      // Calculate abs( delta )
      double delta = iVal - jVal;
      if (delta < 0) delta = -delta;
      Distances.AddElement( delta );
      current++;
    }
  }
}



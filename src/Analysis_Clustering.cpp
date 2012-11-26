// Analysis_Clustering
#include <cfloat> // DBL_MAX
#include "Analysis_Clustering.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // fileExists, integerToString
#include "ProgressBar.h"
#include "DataSet_integer.h" // For converting cnumvtime
#include "Trajout.h"
#include "Analysis_Rms2d.h"

// CONSTRUCTOR
Analysis_Clustering::Analysis_Clustering() :
  epsilon_(-1.0),
  targetNclusters_(-1),
  sieve_(1),
  cnumvtime_(0),
  nofitrms_(false),
  usedme_(false),
  useMass_(false),
  grace_color_(false),
  load_pair_(true),
  cluster_dataset_(0),
  Linkage_(ClusterList::AVERAGELINK),
  clusterfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  singlerepfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  reptrajfmt_(TrajectoryFile::UNKNOWN_TRAJ)
{ } 

void Analysis_Clustering::Help() {
  mprintf("cluster <crd set> [<mask>] [mass] [clusters <n>] [epsilon <e>] [out <cnumvtime>]\n");
  mprintf("        [ linkage | averagelinkage | complete ] [gracecolor] [noload] [nofit] [dme]\n");
  mprintf("        [summary <summaryfile>] [summaryhalf <halffile>] [info <infofile>]\n");
  mprintf("        [ clusterout <trajfileprefix> [clusterfmt <trajformat>] ]\n");
  mprintf("        [ singlerepout <trajfilename> [singlerepfmt <trajformat>] ]\n");
  mprintf("        [ repout <repprefix> [repfmt <repfmt>] ]\n");
  mprintf("        [data <dsetname>]\n");
  mprintf("\tCluster structures based on RMSD or a given DataFile.\n");
  mprintf("<crd set> can be created with the 'createcrd' command.\n");
}

const char* Analysis_Clustering::PAIRDISTFILE = "CpptrajPairDist";

Analysis::RetType Analysis_Clustering::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                                             TopologyList* PFLin, int debugIn)
{
  debug_ = debugIn;
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: clustering: Specify crd set name.\n");
    Help();
    return Analysis::ERR;
  }
  coords_ = (DataSet_Coords*)datasetlist->FindSetOfType( setname, DataSet::COORDS );
  if (coords_ == 0) {
    mprinterr("Error: clustering: Could not locate COORDS set corresponding to %s\n",
              setname.c_str());
    return Analysis::ERR;
  }
  // Check for dataset to cluster on, otherwise coords will be used
  setname = analyzeArgs.GetStringKey("data");
  if (!setname.empty()) {
    cluster_dataset_ = datasetlist->GetDataSet( setname );
    if (cluster_dataset_ == 0) {
      mprinterr("Error: cluster: dataset %s not found.\n", setname.c_str());
      return Analysis::ERR;
    }
  }
  // Get keywords
  useMass_ = analyzeArgs.hasKey("mass");
  targetNclusters_ = analyzeArgs.getKeyInt("clusters",-1);
  sieve_ = analyzeArgs.getKeyInt("sieve",1);
  epsilon_ = analyzeArgs.getKeyDouble("epsilon",-1.0);
  if (analyzeArgs.hasKey("linkage")) Linkage_ = ClusterList::SINGLELINK;
  if (analyzeArgs.hasKey("averagelinkage")) Linkage_ = ClusterList::AVERAGELINK;
  if (analyzeArgs.hasKey("complete")) Linkage_ = ClusterList::COMPLETELINK;
  cnumvtimefile_ = analyzeArgs.GetStringKey("out");
  clusterinfo_ = analyzeArgs.GetStringKey("info");
  summaryfile_ = analyzeArgs.GetStringKey("summary");
  halffile_ = analyzeArgs.GetStringKey("summaryhalf");
  nofitrms_ = analyzeArgs.hasKey("nofit");
  usedme_ = analyzeArgs.hasKey("dme");
  grace_color_ = analyzeArgs.hasKey("gracecolor");
  load_pair_ = !analyzeArgs.hasKey("noload");
  // Output trajectory stuff
  clusterfile_ = analyzeArgs.GetStringKey("clusterout");
  clusterfmt_ = TrajectoryFile::GetFormatFromString( analyzeArgs.GetStringKey("clusterfmt") ); 
  singlerepfile_ = analyzeArgs.GetStringKey("singlerepout");
  singlerepfmt_ = TrajectoryFile::GetFormatFromString( analyzeArgs.GetStringKey("singlerepfmt") );
  reptrajfile_ = analyzeArgs.GetStringKey("repout");
  reptrajfmt_ = TrajectoryFile::GetFormatFromString( analyzeArgs.GetStringKey("repfmt") );
  // Get the mask string 
  maskexpr_ = analyzeArgs.GetMaskNext();

  // Dataset to store cluster number v time
  cnumvtime_ = datasetlist->AddSet(DataSet::INT, analyzeArgs.GetStringNext(), "Cnum");
  if (cnumvtime_==0) return Analysis::ERR;

  // Determine finish criteria. If nothing specified default to 10 clusters.
  if (targetNclusters_==-1 && epsilon_==-1.0)
    targetNclusters_ = 10;

  mprintf("    CLUSTER: Using coords dataset %s, clustering using", coords_->Legend().c_str());
  if ( cluster_dataset_ == 0 ) {
    if (!maskexpr_.empty())
      mprintf(" RMSD (mask [%s])",maskexpr_.c_str());
    else
      mprintf(" RMSD (all atoms)");
    if (useMass_)
      mprintf(", mass-weighted");
    if (nofitrms_)
      mprintf(", no fitting");
    else
      mprintf(" best fit");
  } else 
    mprintf(" dataset %s", cluster_dataset_->Legend().c_str());
  mprintf("\n");
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

  return Analysis::OK;
}

/** This is where the clustering is actually performed. First the distances
  * between each frame are calculated. Then the clustering routine is called.
  */
// TODO: Need to update save to indicate distance type
// NOTE: Should distances be saved only if load_pair?
Analysis::RetType Analysis_Clustering::Analyze() {
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
  if (cluster_dataset_ != 0)
    pairdist_mode = 2;

  if (pairdist_mode == 2) {  // Get distances from dataset.
    calcDistFromDataset( Distances );
    Distances.SaveFile( PAIRDISTFILE );
  }
  if (pairdist_mode == 1) {  // Get distances from file.
    mprintf(" %s found, loading pairwise distances.\n",PAIRDISTFILE);
    if ( Distances.LoadFile(PAIRDISTFILE, coords_->Size()) ) {
      mprintf("\tLoading pairwise distances failed - regenerating from frames.\n");
      pairdist_mode = 0;
    }
  } 
  if (pairdist_mode == 0) { // Get RMSDs between frames
    if (usedme_)
      Analysis_Rms2d::CalcDME(*coords_, Distances, maskexpr_);
    else
      Analysis_Rms2d::Calc2drms(*coords_, Distances, nofitrms_, useMass_, maskexpr_);
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
  return Analysis::OK;
}

// Analysis_Clustering::Print()
void Analysis_Clustering::Print(DataFileList* DFL) {
  // Add dataset to data file list
  DFL->AddSetToFile(cnumvtimefile_, cnumvtime_);

}

// -----------------------------------------------------------------------------
// Analysis_Clustering::ClusterHierAgglo()
/** Cluster using a hierarchical agglomerative (bottom-up) approach. All frames
  * start in their own cluster. The closest two clusters are merged, and 
  * distances between the newly merged cluster and all remaining clusters are
  * recalculated according to one of the following metrics:
  * - single-linkage  : The minimum distance between frames in clusters are used.
  * - average-linkage : The average distance between frames in clusters are used.
  * - complete-linkage: The maximum distance between frames in clusters are used.
  */
int Analysis_Clustering::ClusterHierAgglo( TriangleMatrix& FrameDistances, 
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

// Analysis_Clustering::CreateCnumvtime()
/** Put cluster number vs frame into dataset.  */
void Analysis_Clustering::CreateCnumvtime( ClusterList &CList ) {
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

// Analysis_Clustering::WriteClusterTraj()
/** Write frames in each cluster to a trajectory file.  */
void Analysis_Clustering::WriteClusterTraj( ClusterList &CList ) {
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
    Topology *clusterparm = (Topology*)&(coords_->Top()); // TODO: fix cast
    if (clusterout->SetupTrajWrite(cfilename, 0, clusterparm, clusterfmt_)) 
    {
      mprinterr("Error: Clustering::WriteClusterTraj: Could not set up %s for write.\n",
                cfilename.c_str());
      delete clusterout;
      return;
    }
    //mprinterr("Cluster %i:\n",CList->CurrentNum());
    // Loop over all frames in cluster
    int framenum = 0;
    Frame clusterframe( coords_->Top().Natom() );
    for (; frame != (*C).endframe(); frame++) {
      //mprinterr("%i,",*frame);
      clusterframe.SetFromCRD( (*coords_)[*frame], coords_->NumBoxCrd() );
      clusterout->WriteFrame(framenum++, clusterparm, clusterframe);
    }
    // Close traj
    clusterout->EndTraj();
    //mprinterr("\n");
    //break;
    delete clusterout;
  }
}

// Analysis_Clustering::WriteSingleRepTraj()
/** Write representative frame of each cluster to a trajectory file.  */
void Analysis_Clustering::WriteSingleRepTraj( ClusterList &CList ) {
  Trajout clusterout;

  // Find centroid of first cluster in order to set up parm
  // NOTE: This is redundant if the Summary routine has already been called.
  ClusterList::cluster_iterator cluster = CList.begincluster();
  int framenum = (*cluster).Centroid();

  // Set up trajectory file. Use parm from first frame of cluster (pot. dangerous)
  Topology *clusterparm = (Topology*)&(coords_->Top()); // TODO: fix cast
  if (clusterout.SetupTrajWrite(singlerepfile_, 0, clusterparm, singlerepfmt_)) 
  {
    mprinterr("Error: Clustering::WriteSingleRepTraj: Could not set up %s for write.\n",
                singlerepfile_.c_str());
     return;
  }
  // Write first cluster rep frame
  int framecounter = 0;
  Frame clusterframe( coords_->Top().Natom() );
  clusterframe.SetFromCRD( (*coords_)[ framenum ], coords_->NumBoxCrd() );
  clusterout.WriteFrame(framecounter++, clusterparm, clusterframe);

  ++cluster;
  for (; cluster != CList.endcluster(); cluster++) {
    //mprinterr("Cluster %i: ",CList->CurrentNum());
   framenum = (*cluster).Centroid();
   //mprinterr("%i\n",framenum);
   clusterframe.SetFromCRD( (*coords_)[framenum], coords_->NumBoxCrd() );
   clusterout.WriteFrame(framecounter++, clusterparm, clusterframe);
    //mprinterr("\n");
    //break;
  }
  // Close traj
  clusterout.EndTraj();
}

// Analysis_Clustering::WriteRepTraj()
/** Write representative frame of each cluster to a separate trajectory file,
  * repfile.REPNUM.FMT
  */
void Analysis_Clustering::WriteRepTraj( ClusterList &CList ) {
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
    Topology *clusterparm = (Topology*)&(coords_->Top()); // TODO: Fix cast
    if (clusterout->SetupTrajWrite(cfilename, NULL, clusterparm, reptrajfmt_)) 
    {
      mprinterr("Error: Clustering::WriteRepTraj: Could not set up %s for write.\n",
                reptrajfile_.c_str());
       delete clusterout;
       return;
    }
    // Write cluster rep frame
    Frame clusterframe( coords_->Top().Natom() );
    clusterframe.SetFromCRD( (*coords_)[framenum], coords_->NumBoxCrd() );
    clusterout->WriteFrame(framenum, clusterparm, clusterframe);
    // Close traj
    clusterout->EndTraj();
    delete clusterout;
  }
}

// Analysis_Clustering::calcDistFromDataset()
void Analysis_Clustering::calcDistFromDataset( TriangleMatrix &Distances ) {
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



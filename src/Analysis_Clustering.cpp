// Analysis_Clustering
#include <cfloat> // DBL_MAX
#include "Analysis_Clustering.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // fileExists, integerToString
#include "DataSet_integer.h" // For converting cnumvtime
#include "Trajout.h"

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
  if (sieve_ < 1) {
    mprinterr("Error: 'sieve <#>' must be >= 1 (%i)\n", sieve_);
    return Analysis::ERR;
  }
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
  mprintf("\n\t");
  if (targetNclusters_ != -1)
    mprintf("Looking for %i clusters, ",targetNclusters_);
  if (epsilon_ != -1.0)
    mprintf("Epsilon is %.3f,",epsilon_);
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
  ClusterList CList;
  CList.SetDebug(debug_);

  mprintf("    CLUSTER:");
  // Default: USE_FRAMES  - Calculate pair distances from frames.
  //          USE_FILE    - If PAIRDISTFILE exists, load pair distances from there.
  // Calculated distances will be saved if not loaded from file.
  ClusterList::DistModeType pairdist_mode = ClusterList::USE_FRAMES; 
  if (load_pair_ && fileExists(PAIRDISTFILE))
    pairdist_mode = ClusterList::USE_FILE;
  // If no dataset specified, use COORDS
  if (cluster_dataset_ == 0)
     cluster_dataset_ = (DataSet*) coords_;
  // Calculate distances between frames
  if (CList.CalcFrameDistances( PAIRDISTFILE, cluster_dataset_, pairdist_mode,
                                usedme_, nofitrms_, useMass_, maskexpr_, sieve_ ))
    return Analysis::ERR;
  // Cluster
  CList.ClusterHierAgglo( epsilon_, targetNclusters_, Linkage_);
  // Sort clusters and renumber; also finds centroids for printing
  // representative frames.
  CList.Renumber();
  // If sieving, add remaining frames
  if (sieve_ > 1)
    CList.AddSievedFrames();

  // DEBUG
  if (debug_ > 0) {
    mprintf("\nFINAL CLUSTERS:\n");
    CList.PrintClusters();
    mprintf("\nREPRESENTATIVE FRAMES:\n");
    CList.PrintRepFrames();
  }

  // TEST - print DBI
  mprintf("\tDBI = %f\n", CList.ComputeDBI());

  // Print ptraj-like cluster info
  if (!clusterinfo_.empty())
    CList.PrintClustersToFile(clusterinfo_, coords_->Size());

  // Print a summary of clusters
  if (!summaryfile_.empty())
    CList.Summary(summaryfile_, coords_->Size());

  // Print a summary comparing first half to second half of data for clusters
  if (!halffile_.empty())
    CList.Summary_Half(halffile_, coords_->Size());

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
// Analysis_Clustering::CreateCnumvtime()
/** Put cluster number vs frame into dataset.  */
void Analysis_Clustering::CreateCnumvtime( ClusterList &CList ) {
  // FIXME:
  // Cast generic DataSet for cnumvtime back to integer dataset to 
  // access specific integer dataset functions for resizing and []
  // operator. Should this eventually be generic to all atomic DataSets? 
  DataSet_integer* cnum_temp = (DataSet_integer*)cnumvtime_;
  cnum_temp->Resize( coords_->Size() );

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
      coords_->GetFrame( *frame, clusterframe );
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
  // Set up trajectory file. Use parm from COORDS DataSet. 
  Topology *clusterparm = (Topology*)&(coords_->Top()); // TODO: fix cast
  if (clusterout.SetupTrajWrite(singlerepfile_, 0, clusterparm, singlerepfmt_)) 
  {
    mprinterr("Error: Clustering::WriteSingleRepTraj: Could not set up %s for write.\n",
                singlerepfile_.c_str());
     return;
  }
  // Set up frame to hold cluster rep coords. 
  Frame clusterframe( coords_->Top().Natom() );
  int framecounter = 0;
  // Write rep frames from all clusters.
  for (ClusterList::cluster_iterator cluster = CList.begincluster(); 
                                     cluster != CList.endcluster(); ++cluster) 
  {
   coords_->GetFrame( (*cluster).CentroidFrame(), clusterframe );
   clusterout.WriteFrame(framecounter++, clusterparm, clusterframe);
  }
  // Close traj
  clusterout.EndTraj();
}

// Analysis_Clustering::WriteRepTraj()
/** Write representative frame of each cluster to a separate trajectory file,
  * repfile.REPNUM.FMT
  */
void Analysis_Clustering::WriteRepTraj( ClusterList &CList ) {
  // Get extension for representative frame format 
  std::string tmpExt = TrajectoryFile::GetExtensionForType(reptrajfmt_);
  // Use Topology from COORDS DataSet to set up input frame
  Topology *clusterparm = (Topology*)&(coords_->Top()); // TODO: Fix cast
  Frame clusterframe( clusterparm->Natom() );
  // Loop over all clusters
  for (ClusterList::cluster_iterator C = CList.begincluster();
                                     C != CList.endcluster(); ++C)
  {
    Trajout* clusterout = new Trajout();
    // Get centroid frame # 
    int framenum = (*C).CentroidFrame();
    // Create filename based on frame #
    std::string cfilename = reptrajfile_ + "." + integerToString(framenum+1) + tmpExt;
    // Set up trajectory file. 
    if (clusterout->SetupTrajWrite(cfilename, 0, clusterparm, reptrajfmt_)) 
    {
      mprinterr("Error: Clustering::WriteRepTraj: Could not set up %s for write.\n",
                cfilename.c_str());
       delete clusterout;
       return;
    }
    // Write cluster rep frame
    coords_->GetFrame( framenum, clusterframe );
    clusterout->WriteFrame(framenum, clusterparm, clusterframe);
    // Close traj
    clusterout->EndTraj();
    delete clusterout;
  }
}

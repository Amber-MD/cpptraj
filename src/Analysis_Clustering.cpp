// Analysis_Clustering
#include "Analysis_Clustering.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // fileExists, integerToString
#include "DataSet_integer.h" // For converting cnumvtime
#include "Trajout.h"
// Clustering Algorithms
#include "Cluster_HierAgglo.h"
#include "Cluster_DBSCAN.h"

// CONSTRUCTOR
Analysis_Clustering::Analysis_Clustering() :
  masterDSL_(0),
  coords_(0),
  CList_(0),
  sieve_(1),
  splitFrame_(-1),
  cnumvtime_(0),
  cpopvtimefile_(0),
  nofitrms_(false),
  usedme_(false),
  useMass_(false),
  grace_color_(false),
  norm_pop_(false),
  load_pair_(false),
  clusterfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  singlerepfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  reptrajfmt_(TrajectoryFile::UNKNOWN_TRAJ)
{ } 

// DESTRUCTOR
Analysis_Clustering::~Analysis_Clustering() {
  if (CList_ != 0) delete CList_;
}

void Analysis_Clustering::Help() {
  mprintf("\t[crdset <crd set>]\n");
  mprintf("  Algorithms:\n");
  Cluster_HierAgglo::Help();
  Cluster_DBSCAN::Help();
  mprintf("  Distance options:\n");
  mprintf("\t{[[rms] [<mask>] [mass] [nofit]] | [dme [<mask>]] | [data <dset0>[,<dset1>,...]]}\n");
  mprintf("\t[sieve <#>] [loadpairdist] [savepairdist] [pairdist <file>]\n");
  mprintf("  Output options:\n");
  mprintf("\t[out <cnumvtime>] [gracecolor] [summary <summaryfile>] [info <infofile>]\n");
  mprintf("\t[summaryhalf <halffile>] [splitframe <frame>] [cpopvtime <file> [normpop]]\n");
  mprintf("  Coordinate output options:\n");
  mprintf("\t[ clusterout <trajfileprefix> [clusterfmt <trajformat>] ]\n");
  mprintf("\t[ singlerepout <trajfilename> [singlerepfmt <trajformat>] ]\n");
  mprintf("\t[ repout <repprefix> [repfmt <repfmt>] ]\n");
  mprintf("\tCluster structures based on coordinates (RMSD/DME) or given data set(s).\n");
  mprintf("\t<crd set> can be created with the 'createcrd' command.\n");
}

const char* Analysis_Clustering::PAIRDISTFILE = "CpptrajPairDist";

Analysis::RetType Analysis_Clustering::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  debug_ = debugIn;
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringKey("crdset");
  coords_ = (DataSet_Coords*)datasetlist->FindCoordsSet( setname );
  if (coords_ == 0) {
    mprinterr("Error: clustering: Could not locate COORDS set corresponding to %s\n",
              setname.c_str());
    Help();
    return Analysis::ERR;
  }
  // Check for DataSet(s) to cluster on, otherwise coords will be used
  cluster_dataset_.clear();
  setname = analyzeArgs.GetStringKey("data");
  if (!setname.empty()) {
    ArgList dsnames(setname, ",");
    for (ArgList::const_iterator name = dsnames.begin(); name != dsnames.end(); ++name) {
      DataSet* ds = datasetlist->GetDataSet( *name );
      if (ds == 0) {
        mprinterr("Error: cluster: dataset %s not found.\n", (*name).c_str());
        return Analysis::ERR;
      }
      cluster_dataset_.push_back( ds );
    }
  }
  // Get clustering algorithm
  if (CList_ != 0) delete CList_;
  CList_ = 0;
  if (analyzeArgs.hasKey("hieragglo"))   CList_ = new Cluster_HierAgglo(); 
  else if (analyzeArgs.hasKey("dbscan")) CList_ = new Cluster_DBSCAN();
  else {
    mprintf("Warning: No clustering algorithm specified; defaulting to 'hieragglo'\n");
    CList_ = new Cluster_HierAgglo();
  }
  if (CList_ == 0) return Analysis::ERR;
  CList_->SetDebug(debug_);
  // Get algorithm-specific keywords
  if (CList_->SetupCluster( analyzeArgs )) return Analysis::ERR; 
  // Get keywords
  useMass_ = analyzeArgs.hasKey("mass");
  usedme_ = analyzeArgs.hasKey("dme");
  sieve_ = analyzeArgs.getKeyInt("sieve",1);
  if (sieve_ < 1) {
    mprinterr("Error: 'sieve <#>' must be >= 1 (%i)\n", sieve_);
    return Analysis::ERR;
  }
  splitFrame_ = analyzeArgs.getKeyInt("splitframe", -1);
  DataFile* cnumvtimefile = DFLin->AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  cpopvtimefile_ = DFLin->AddDataFile(analyzeArgs.GetStringKey("cpopvtime"), analyzeArgs);
  clusterinfo_ = analyzeArgs.GetStringKey("info");
  summaryfile_ = analyzeArgs.GetStringKey("summary");
  halffile_ = analyzeArgs.GetStringKey("summaryhalf");
  nofitrms_ = analyzeArgs.hasKey("nofit");
  grace_color_ = analyzeArgs.hasKey("gracecolor");
  if (cpopvtimefile_ != 0) {
    if (grace_color_) {
      mprintf("Warning: 'gracecolor' not compatible with 'cpopvtime' - disabling 'gracecolor'\n");
      grace_color_ = false;
    }
    norm_pop_ = analyzeArgs.hasKey("normpop");
  }
  // Options for loading/saving pairwise distance file
  load_pair_ = analyzeArgs.hasKey("loadpairdist");
  bool save_pair = analyzeArgs.hasKey("savepairdist");
  pairdistfile_ = analyzeArgs.GetStringKey("pairdist");
  if ( (load_pair_ || save_pair) && pairdistfile_.empty() )
    pairdistfile_.assign(PAIRDISTFILE);
  else if (!pairdistfile_.empty())
    load_pair_ = true;
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
  if (cnumvtimefile != 0) cnumvtimefile->AddSet( cnumvtime_ ); 
  // Save master DSL for Cpopvtime
  masterDSL_ = datasetlist;

  mprintf("    CLUSTER: Using coords dataset %s, clustering using", coords_->Legend().c_str());
  if ( cluster_dataset_.empty() ) {
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
  } else {
    if (cluster_dataset_.size() == 1)
      mprintf(" dataset %s", cluster_dataset_[0]->Legend().c_str());
    else
      mprintf(" %u datasets.", cluster_dataset_.size());
  }
  mprintf("\n");
  CList_->ClusteringInfo();
  if (sieve_ > 1)
    mprintf("\tInitial clustering sieve value is %i frames.\n", sieve_);
  if (cnumvtimefile != 0)
    mprintf("\tCluster # vs time will be written to %s\n", cnumvtimefile->DataFilename().base());
  if (cpopvtimefile_ != 0) {
    mprintf("\tCluster pop vs time will be written to %s", cpopvtimefile_->DataFilename().base());
    if (norm_pop_) mprintf(" (normalized)");
    mprintf("\n");
  }
  if (grace_color_)
    mprintf("\tGrace color instead of cluster number (1-15) will be saved.\n");
  if (load_pair_)
    mprintf("\tPreviously calcd pair distances %s will be used if found.\n",
            pairdistfile_.c_str());
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

  return Analysis::OK;
}

/** This is where the clustering is actually performed. First the distances
  * between each frame are calculated. Then the clustering routine is called.
  */
// TODO: Need to update save to indicate distance type
// NOTE: Should distances be saved only if load_pair?
Analysis::RetType Analysis_Clustering::Analyze() {
  mprintf("\tStarting clustering.\n");
  // Default: USE_FRAMES  - Calculate pair distances from frames.
  //          USE_FILE    - If pairdistfile exists, load pair distances from there.
  // Calculated distances will be saved if not loaded from file.
  ClusterList::DistModeType pairdist_mode = ClusterList::USE_FRAMES; 
  if (load_pair_ && fileExists(pairdistfile_.c_str()))
    pairdist_mode = ClusterList::USE_FILE;
  // If no dataset specified, use COORDS
  if (cluster_dataset_.empty())
     cluster_dataset_.push_back( (DataSet*)coords_ );
  // Calculate distances between frames
  if (CList_->CalcFrameDistances( pairdistfile_, cluster_dataset_, pairdist_mode,
                                  usedme_, nofitrms_, useMass_, maskexpr_, sieve_ ))
    return Analysis::ERR;
  // Cluster
  CList_->Cluster();
  // Sort clusters and renumber; also finds centroids for printing
  // representative frames.
  CList_->Renumber();
  // If sieving, add remaining frames
  if (sieve_ > 1)
    CList_->AddSievedFrames();

  // DEBUG
  if (debug_ > 0) {
    mprintf("\nFINAL CLUSTERS:\n");
    CList_->PrintClusters();
  }

  // Print ptraj-like cluster info
  if (!clusterinfo_.empty())
    CList_->PrintClustersToFile(clusterinfo_, coords_->Size());

  // Print a summary of clusters
  if (!summaryfile_.empty())
    CList_->Summary(summaryfile_, coords_->Size());

  // Print a summary comparing first half to second half of data for clusters
  if (!halffile_.empty())
    CList_->Summary_Half(halffile_, coords_->Size(), splitFrame_);

  // Create cluster v time data from clusters.
  CreateCnumvtime( *CList_ );

  // Create cluster pop v time plots
  if (cpopvtimefile_ != 0)
    CreateCpopvtime( *CList_ );

  // Write clusters to trajectories
  if (!clusterfile_.empty())
    WriteClusterTraj( *CList_ ); 

  // Write all representative frames to a single traj
  if (!singlerepfile_.empty())
    WriteSingleRepTraj( *CList_ );

  // Write all representative frames to separate trajs
  if (!reptrajfile_.empty())
    WriteRepTraj( *CList_ );
  return Analysis::OK;
}

// -----------------------------------------------------------------------------
// Analysis_Clustering::CreateCnumvtime()
/** Put cluster number vs frame into dataset.  */
void Analysis_Clustering::CreateCnumvtime( ClusterList const& CList ) {
  // FIXME:
  // Cast generic DataSet for cnumvtime back to integer dataset to 
  // access specific integer dataset functions for resizing and []
  // operator. Should this eventually be generic to all atomic DataSets? 
  DataSet_integer* cnum_temp = (DataSet_integer*)cnumvtime_;
  cnum_temp->Resize( coords_->Size() );
  // Make all clusters start at -1. This way cluster algorithms that
  // have noise points (i.e. no cluster assigned) will be distinguished.
  if (!grace_color_)
    for (int i = 0; i < cnum_temp->Size(); ++i)
      (*cnum_temp)[i] = -1;

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

// Analysis_Clustering::CreateCpopvtime()
// NOTE: Should not be called if cpopvtimefile is NULL
void Analysis_Clustering::CreateCpopvtime( ClusterList const& CList ) {
  std::vector<int> Pop(CList.Nclusters(), 0);
  // Set up normalization
  std::vector<double> Norm(CList.Nclusters(), 1.0);
  if (norm_pop_) {
    int cnum = 0;
    for (ClusterList::cluster_iterator C = CList.begincluster(); 
                                       C != CList.endcluster(); ++C)
      Norm[cnum++] = 1.0 / (double)((*C).Nframes());
  }
  std::vector<DataSet*> DSL;
  for (int cnum = 0; cnum < CList.Nclusters(); ++cnum) { 
    DSL.push_back(masterDSL_->AddSetIdxAspect( DataSet::FLOAT, cnumvtime_->Name(), 
                                               cnum, "Pop" ));
    if (DSL.back() == 0) {
      mprinterr("Error: Could not allocate cluster pop v time DataSet.\n");
      return;
    }
    cpopvtimefile_->AddSet( DSL.back() );
  }
  // Assumes cnumvtime has been calcd and not gracecolor!
  // TODO: Normalization
  DataSet_integer* cnum_temp = (DataSet_integer*)cnumvtime_;
  for (int frame = 0; frame < coords_->Size(); ++frame) {
    int cluster_num = (*cnum_temp)[frame];
    // Noise points are -1
    if (cluster_num > -1)
      Pop[cluster_num]++;
    for (int cnum = 0; cnum < CList.Nclusters(); ++cnum) {
      float f = ((double)Pop[cnum] * Norm[cnum]);
      DSL[cnum]->Add(frame, &f);
    }
  }
}

// Analysis_Clustering::WriteClusterTraj()
/** Write frames in each cluster to a trajectory file.  */
void Analysis_Clustering::WriteClusterTraj( ClusterList const& CList ) {
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
void Analysis_Clustering::WriteSingleRepTraj( ClusterList const& CList ) {
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
void Analysis_Clustering::WriteRepTraj( ClusterList const& CList ) {
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

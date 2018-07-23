// Analysis_Clustering
#include "Analysis_Clustering.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // fileExists, integerToString
#include "DataSet_integer.h" // For converting cnumvtime
#include "DataSet_float.h"
#include "Trajout_Single.h"
#include "Timer.h"
// Clustering Algorithms
#include "Cluster_HierAgglo.h"
#include "Cluster_DBSCAN.h"
#include "Cluster_Kmeans.h"
#include "Cluster_ReadInfo.h"
#include "Cluster_DPeaks.h"

// CONSTRUCTOR
Analysis_Clustering::Analysis_Clustering() :
  masterDSL_(0),
  coords_(0),
  CList_(0),
  sieve_(1),
  sieveSeed_(-1),
  windowSize_(0),
  drawGraph_(0),
  draw_maxit_(0),
  nRepsToSave_(1),
  draw_tol_(0.0),
  refCut_(1.0),
  cnumvtime_(0),
  clustersVtime_(0),
  pw_dist_(0),
  cpopvtimefile_(0),
  pwd_file_(0),
  nofitrms_(false),
  metric_(ClusterList::RMS),
  useMass_(false),
  grace_color_(false),
  norm_pop_(NONE),
  bestRep_(CUMULATIVE),
  calc_lifetimes_(false),
  writeRepFrameNum_(false),
  includeSieveInCalc_(false),
  clusterfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  singlerepfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  reptrajfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  debug_(0)
{ } 

const TrajectoryFile::TrajFormatType Analysis_Clustering::DEF_TRAJ_FMT_ = TrajectoryFile::AMBERTRAJ;

// DESTRUCTOR
Analysis_Clustering::~Analysis_Clustering() {
  if (CList_ != 0) delete CList_;
}

void Analysis_Clustering::Help() const {
  mprintf("\t[crdset <crd set> | nocoords]\n");
  mprintf("  Algorithms:\n");
  Cluster_HierAgglo::Help();
  Cluster_DBSCAN::Help();
  Cluster_DPeaks::Help();
  Cluster_Kmeans::Help();
  Cluster_ReadInfo::Help();
  mprintf("  Distance metric options: {rms | srmsd | dme | data}\n"
          "\t{ [[rms | srmsd] [<mask>] [mass] [nofit]] | [dme [<mask>]] |\n"
          "\t   [data <dset0>[,<dset1>,...]] }\n"
          "\t[sieve <#> [random [sieveseed <#>]]] [loadpairdist] [savepairdist] [pairdist <name>]\n"
          "\t[pairwisecache {mem | disk | none}] [includesieveincalc]\n"
          "  Output options:\n"
          "\t[out <cnumvtime>] [gracecolor] [summary <summaryfile>] [info <infofile>]\n"
          "\t[summarysplit <splitfile>] [splitframe <comma-separated frame list>]\n"
          "\t[bestrep {cumulative|centroid|cumulative_nosieve}] [savenreps <#>]\n"
          "\t[clustersvtime <filename> cvtwindow <window size>]\n"
          "\t[cpopvtime <file> [normpop | normframe]] [lifetime]\n"
          "\t[sil <silhouette file prefix>] [assignrefs [refcut <rms>] [refmask <mask>]]\n"
          "  Coordinate output options:\n"
          "\t[ clusterout <trajfileprefix> [clusterfmt <trajformat>] ]\n"
          "\t[ singlerepout <trajfilename> [singlerepfmt <trajformat>] ]\n"
          "\t[ repout <repprefix> [repfmt <repfmt>] [repframe] ]\n"
          "\t[ avgout <avgprefix> [avgfmt <avgfmt>] ]\n"
          "  Experimental options:\n"
          "\t[[drawgraph | drawgraph3d] [draw_tol <tolerance>] [draw_maxit <iterations]]\n"
          "  Cluster structures based on coordinates (RMSD/DME) or given data set(s).\n"
          "  <crd set> can be created with the 'createcrd' command.\n");
          /// pytraj can turn off cluster info by specifying 'noinfo' keyword
}

const char* Analysis_Clustering::PAIRDISTFILE_ = "CpptrajPairDist";
DataFile::DataFormatType Analysis_Clustering::PAIRDISTTYPE_ =
# ifdef BINTRAJ
  DataFile::NCCMATRIX;
# else
  DataFile::CMATRIX;
# endif

// Analysis_Clustering::GetClusterTrajArgs()
void Analysis_Clustering::GetClusterTrajArgs(ArgList& argIn,
                                             const char* trajKey, const char* fmtKey,
                                             std::string& trajName,
                                             TrajectoryFile::TrajFormatType& fmt) const
{
  trajName = argIn.GetStringKey( trajKey );
  fmt = TrajectoryFile::WriteFormatFromString( argIn.GetStringKey(fmtKey), fmt );
  // If file name specified but not format, try to guess from name
  if (!trajName.empty() && fmt == TrajectoryFile::UNKNOWN_TRAJ)
    fmt = TrajectoryFile::WriteFormatFromFname( trajName, DEF_TRAJ_FMT_ );
}

// Analysis_Clustering::Setup()
Analysis::RetType Analysis_Clustering::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  debug_ = debugIn;
  if (analyzeArgs.hasKey("nocoords"))
    coords_ = 0;
  else {
    // Attempt to get coords dataset from datasetlist
    std::string setname = analyzeArgs.GetStringKey("crdset");
    coords_ = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
    if (coords_ == 0) {
      mprinterr("Error: Could not locate COORDS set corresponding to %s\n",
                setname.c_str());
      return Analysis::ERR;
    }
  }
  // Check for DataSet(s) to cluster on, otherwise coords will be used
  cluster_dataset_.clear();
  std::string dataSetname = analyzeArgs.GetStringKey("data");
  metric_ = ClusterList::RMS;
  if (!dataSetname.empty()) {
    ArgList dsnames(dataSetname, ",");
    DataSetList inputDsets;
    for (ArgList::const_iterator name = dsnames.begin(); name != dsnames.end(); ++name) {
      DataSetList tempDSL = setup.DSL().GetMultipleSets( *name );
      if (tempDSL.empty()) {
        mprinterr("Error: %s did not correspond to any data sets.\n", dataSetname.c_str());
        return Analysis::ERR;
      }
      inputDsets += tempDSL;
    }
    for (DataSetList::const_iterator ds = inputDsets.begin(); ds != inputDsets.end(); ++ds) {
      // Clustering only allowed on 1D data sets.
      if ( (*ds)->Ndim() != 1 ) {
        mprinterr("Error: Clustering only allowed on 1D data sets, %s is %zuD.\n",
                  (*ds)->legend(), (*ds)->Ndim());
        return Analysis::ERR;
      }
      cluster_dataset_.push_back( *ds );
    }
    metric_ = ClusterList::DATA;
  } else {
    int usedme = (int)analyzeArgs.hasKey("dme");
    int userms = (int)analyzeArgs.hasKey("rms");
    int usesrms = (int)analyzeArgs.hasKey("srmsd");
    if (usedme + userms + usesrms > 1) {
      mprinterr("Error: Specify either 'dme', 'rms', or 'srmsd'.\n");
      return Analysis::ERR;
    }
    if      (usedme)  metric_ = ClusterList::DME;
    else if (userms)  metric_ = ClusterList::RMS;
    else if (usesrms) metric_ = ClusterList::SRMSD;
  }
  // Get all loaded reference structures
  if (analyzeArgs.hasKey("assignrefs")) {
    refs_ = setup.DSL().GetSetsOfType("*", DataSet::REF_FRAME);
    if (refs_.empty()) {
      mprinterr("Error: 'assignrefs' specified but no references loaded.\n");
      return Analysis::ERR;
    }
    refCut_ = analyzeArgs.getKeyDouble("refcut", 1.0);
    refmaskexpr_ = analyzeArgs.GetStringKey("refmask");
  }
  // Get clustering algorithm
  if (CList_ != 0) delete CList_;
  CList_ = 0;
  if (analyzeArgs.hasKey("hieragglo"))   CList_ = new Cluster_HierAgglo(); 
  else if (analyzeArgs.hasKey("dbscan")) CList_ = new Cluster_DBSCAN();
  else if (analyzeArgs.hasKey("dpeaks")) CList_ = new Cluster_DPeaks();
  else if (analyzeArgs.hasKey("kmeans") ||
           analyzeArgs.hasKey("means" )) CList_ = new Cluster_Kmeans();
  else if (analyzeArgs.hasKey("readinfo") ||
           analyzeArgs.hasKey("readtxt")) CList_ = new Cluster_ReadInfo(); 
  else {
    mprintf("Warning: No clustering algorithm specified; defaulting to 'hieragglo'\n");
    CList_ = new Cluster_HierAgglo();
  }
  if (CList_ == 0) return Analysis::ERR;
  CList_->SetDebug(debug_);
  // Get algorithm-specific keywords
  if (CList_->SetupCluster( analyzeArgs )) return Analysis::ERR; 
  // Get keywords
  includeSieveInCalc_ = analyzeArgs.hasKey("includesieveincalc");
  if (includeSieveInCalc_)
    mprintf("Warning: 'includesieveincalc' may be very slow.\n");
  useMass_ = analyzeArgs.hasKey("mass");
  sieveSeed_ = analyzeArgs.getKeyInt("sieveseed", -1);
  sieve_ = analyzeArgs.getKeyInt("sieve", 1);
  if (sieve_ < 1) {
    mprinterr("Error: 'sieve <#>' must be >= 1 (%i)\n", sieve_);
    return Analysis::ERR;
  }
  if (analyzeArgs.hasKey("random") && sieve_ > 1)
    sieve_ = -sieve_; // negative # indicates random sieve
  halffile_ = analyzeArgs.GetStringKey("summarysplit");
  if (halffile_.empty()) // For backwards compat.
    halffile_ = analyzeArgs.GetStringKey("summaryhalf");
  if (!halffile_.empty()) {
    ArgList splits( analyzeArgs.GetStringKey("splitframe"), "," );
    if (!splits.empty()) {
      splitFrames_.clear();
      int sf = splits.getNextInteger(-1); // User frame #s start at 1
      while (sf > 0) {
        splitFrames_.push_back( sf );
        sf = splits.getNextInteger(-1);
      }
      if ((int)splitFrames_.size() < splits.Nargs()) {
        mprinterr("Error: Invalid split frame arguments.\n");
        splits.CheckForMoreArgs();
        return Analysis::ERR;
      }
    }
  }
  std::string bestRepStr = analyzeArgs.GetStringKey("bestrep");
  if (bestRepStr.empty()) {
    // For sieving, cumulative can get very expensive. Default to centroid.
    if (sieve_ != 1)
      bestRep_ = CENTROID;
    else
      bestRep_ = CUMULATIVE;
  } else {
    if (bestRepStr == "cumulative")
      bestRep_ = CUMULATIVE;
    else if (bestRepStr == "centroid")
      bestRep_ = CENTROID;
    else if (bestRepStr == "cumulative_nosieve")
      bestRep_ = CUMULATIVE_NOSIEVE;
    else {
      mprinterr("Error: Invalid 'bestRep' option (%s)\n", bestRepStr.c_str());
      return Analysis::ERR;
    }
  }
  nRepsToSave_ = analyzeArgs.getKeyInt("savenreps", 1);
  if (nRepsToSave_ < 1) {
    mprinterr("Error: 'savenreps' must be > 0\n");
    return Analysis::ERR;
  }
  if (analyzeArgs.hasKey("drawgraph"))
    drawGraph_ = 1;
  else if (analyzeArgs.hasKey("drawgraph3d"))
    drawGraph_ = 2;
  else
    drawGraph_ = 0;
  draw_maxit_ = analyzeArgs.getKeyInt("draw_maxit", 1000);
  draw_tol_ = analyzeArgs.getKeyDouble("draw_tol", 1.0E-5);
  
  DataFile* cnumvtimefile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  DataFile* clustersvtimefile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("clustersvtime"),
                                                   analyzeArgs);
  windowSize_ = analyzeArgs.getKeyInt("cvtwindow", 0);
  cpopvtimefile_ = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("cpopvtime"), analyzeArgs);
  clusterinfo_ = analyzeArgs.GetStringKey("info");
  summaryfile_ = analyzeArgs.GetStringKey("summary");
  nofitrms_ = analyzeArgs.hasKey("nofit");
  grace_color_ = analyzeArgs.hasKey("gracecolor");
  calc_lifetimes_ = analyzeArgs.hasKey("lifetime");
  if (cpopvtimefile_ != 0) {
    if (analyzeArgs.hasKey("normpop"))
      norm_pop_ = CLUSTERPOP;
    else if (analyzeArgs.hasKey("normframe"))
      norm_pop_ = FRAME;
    else
      norm_pop_ = NONE;
  }
  sil_file_ = analyzeArgs.GetStringKey("sil");
  // ---------------------------------------------
  // Options for loading/saving pairwise distance file
  DataSet::DataType pw_type = DataSet::CMATRIX;
  std::string pw_typeString = analyzeArgs.GetStringKey("pairwisecache");
  if (!pw_typeString.empty()) {
    if (pw_typeString == "mem")
      pw_type = DataSet::CMATRIX;
    else if (pw_typeString == "disk")
      pw_type = DataSet::CMATRIX_DISK;
    else if (pw_typeString == "none")
      pw_type = DataSet::CMATRIX_NOMEM;
    else {
      mprinterr("Error: Unrecognized option for 'pairwisecache' ('%s')\n", pw_typeString.c_str());
      return Analysis::ERR;
    }
  }
  std::string pairdistname = analyzeArgs.GetStringKey("pairdist");
  DataFile::DataFormatType pairdisttype = DataFile::UNKNOWN_DATA;
  bool load_pair = analyzeArgs.hasKey("loadpairdist");
  bool save_pair = analyzeArgs.hasKey("savepairdist");
  pw_dist_ = 0;
  if (load_pair) {
    // If 'loadpairdist' specified, assume we want to load from file.
    if (pairdistname.empty()) {
      pairdistname = PAIRDISTFILE_;
      pairdisttype = PAIRDISTTYPE_;
    }
    if (File::Exists( pairdistname )) {
      DataFile dfIn;
      if (dfIn.ReadDataIn( pairdistname, ArgList(), setup.DSL() )) return Analysis::ERR;
      pw_dist_ = setup.DSL().GetDataSet( pairdistname );
      if (pw_dist_ == 0) return Analysis::ERR;
    } else
      pairdisttype = PAIRDISTTYPE_;
  }
  if (pw_dist_ == 0 && !pairdistname.empty()) {
    // Just 'pairdist' specified or loadpairdist specified and file not found.
    // Look for Cmatrix data set. 
    pw_dist_ = setup.DSL().FindSetOfType( pairdistname, DataSet::CMATRIX );
    //if (pw_dist_ == 0) { // TODO: Convert general matrix to cluster matrix
    //  mprinterr("Error: Cluster matrix with name '%s' not found.\n");
    //  return Analysis::ERR;
    //}
    if (pw_dist_ == 0 && load_pair) {
    // If the file (or dataset) does not yet exist we will assume we want to save.
      mprintf("Warning: 'loadpairdist' specified but '%s' not found; will save distances.\n",
              pairdistname.c_str());
      save_pair = true;
    }
  }
  // Create file for saving pairwise distances
  pwd_file_ = 0;
  if (save_pair) {
    if (pairdistname.empty()) {
      pairdistname = PAIRDISTFILE_;
      pairdisttype = PAIRDISTTYPE_;
    }
    pwd_file_ = setup.DFL().AddDataFile( pairdistname, pairdisttype, ArgList() );
  }
  // ---------------------------------------------
  // Output trajectory stuff
  writeRepFrameNum_ = analyzeArgs.hasKey("repframe");
  GetClusterTrajArgs(analyzeArgs, "clusterout",   "clusterfmt",   clusterfile_,   clusterfmt_);
  GetClusterTrajArgs(analyzeArgs, "singlerepout", "singlerepfmt", singlerepfile_, singlerepfmt_);
  GetClusterTrajArgs(analyzeArgs, "repout",       "repfmt",       reptrajfile_,   reptrajfmt_);
  GetClusterTrajArgs(analyzeArgs, "avgout",       "avgfmt",       avgfile_,       avgfmt_);

  // Get the mask string 
  maskexpr_ = analyzeArgs.GetMaskNext();
  if (!refs_.empty() && refmaskexpr_.empty()) {
    refmaskexpr_ = maskexpr_;
    if (refmaskexpr_.empty()) {
      refmaskexpr_.assign("!@H=");
      mprintf("Warning: 'assignrefs' specified but no 'refmask' given.\n"
              "Warning:   Using default mask expression: '%s'\n", refmaskexpr_.c_str());
    }
  }

  // Output option for cluster info
  suppressInfo_ = analyzeArgs.hasKey("noinfo");

  // Dataset to store cluster number v time
  cnumvtime_ = setup.DSL().AddSet(DataSet::INTEGER, analyzeArgs.GetStringNext(), "Cnum");
  if (cnumvtime_==0) return Analysis::ERR;
  if (cnumvtimefile != 0) cnumvtimefile->AddDataSet( cnumvtime_ );
  // If no pairwise distance matrix yet, allocate one
  if (pw_dist_ == 0) {
    MetaData md;
    if (pairdistname.empty())
      md = MetaData( cnumvtime_->Meta().Name(), "PWD" );
    else
      md = MetaData( pairdistname );
    if (pw_type == DataSet::CMATRIX_DISK)
      md.SetFileName("CpptrajPairwiseCache");
    pw_dist_ = setup.DSL().AddSet(pw_type, md);
    if (pw_dist_ == 0) return Analysis::ERR;
  }

  // DataSet for # clusters seen v time
  if (clustersvtimefile != 0) {
    if (windowSize_ < 2) {
      mprinterr("Error: For # clusters seen vs time, cvtwindow must be specified and > 1\n");
      return Analysis::ERR;
    }
    clustersVtime_ = setup.DSL().AddSet(DataSet::INTEGER, 
                                         MetaData(cnumvtime_->Meta().Name(), "NCVT"));
    if (clustersVtime_ == 0) return Analysis::ERR;
    clustersvtimefile->AddDataSet( clustersVtime_ );
  }
  // Save master DSL for Cpopvtime
  masterDSL_ = setup.DslPtr();

  mprintf("    CLUSTER:");
  if (coords_ != 0) mprintf(" Using coords dataset %s,", coords_->legend());
  mprintf(" clustering using");
  if ( metric_ != ClusterList::DATA ) {
    mprintf(" %s", ClusterList::MetricString( metric_ ));
    if (!maskexpr_.empty())
      mprintf(" (mask [%s])",maskexpr_.c_str());
    else
      mprintf(" (all atoms)");
    if (useMass_)
      mprintf(", mass-weighted");
    if (nofitrms_)
      mprintf(", no fitting");
    else
      mprintf(" best-fit");
  } else {
    if (cluster_dataset_.size() == 1)
      mprintf(" dataset %s", cluster_dataset_[0]->legend());
    else
      mprintf(" %u datasets.", cluster_dataset_.size());
  }
  mprintf("\n");
  CList_->ClusteringInfo();
  if (sieve_ > 1)
    mprintf("\tInitial clustering sieve value is %i frames.\n", sieve_);
  else if (sieve_ < -1) {
    mprintf("\tInitial clustering will be randomly sieved (with value %i)", -sieve_);
    if (sieveSeed_ > 0) mprintf(" using random seed %i", sieveSeed_);
    mprintf(".\n");
  }
  if (sieve_ != 1) {
    if (includeSieveInCalc_)
      mprintf("\tAll frames (including sieved) will be used to calc within-cluster average.\n");
    else
      mprintf("\tOnly non-sieved frames will be used to calc within-cluster average.\n");
  }
  if (cnumvtimefile != 0)
    mprintf("\tCluster # vs time will be written to %s\n", cnumvtimefile->DataFilename().base());
  if (clustersvtimefile != 0)
    mprintf("\t# clusters seen vs time will be written to %s\n",
            clustersvtimefile->DataFilename().base());
  if (cpopvtimefile_ != 0) {
    mprintf("\tCluster pop vs time will be written to %s", cpopvtimefile_->DataFilename().base());
    if (norm_pop_==CLUSTERPOP) mprintf(" (normalized by cluster size)");
    else if (norm_pop_==FRAME) mprintf(" (normalized by frame)");
    mprintf("\n");
  }
  if (grace_color_)
    mprintf("\tGrace color instead of cluster number (1-15) will be saved.\n");
  if (calc_lifetimes_)
    mprintf("\tCluster lifetime data sets will be calculated.\n");
  mprintf("\tPairwise distance data set is '%s'\n", pw_dist_->legend());
  if (pw_dist_->Type() == DataSet::CMATRIX_NOMEM)
    mprintf("\tPairwise distances will not be cached (will slow clustering calcs)\n");
  else if (pw_dist_->Type() == DataSet::CMATRIX_DISK)
    mprintf("\tPairwise distances will be cached to disk (will slow clustering calcs)\n");
  if (pwd_file_ != 0)
    mprintf("\tSaving pair-wise distances to '%s'\n", pwd_file_->DataFilename().full());
  if (!clusterinfo_.empty())
    mprintf("\tCluster information will be written to %s\n",clusterinfo_.c_str());
  if (!summaryfile_.empty())
    mprintf("\tSummary of cluster results will be written to %s\n",summaryfile_.c_str());
  if (!sil_file_.empty()) {
    mprintf("\tFrame silhouettes will be written to %s.frame.dat, cluster silhouettes\n"
            "\t  will be written to %s.cluster.dat\n", sil_file_.c_str(), sil_file_.c_str());
    if (sieve_ != 1) {
      if (includeSieveInCalc_)
        mprintf("\tSilhouette calculation will use all frames.\n");
      else
        mprintf("\tSilhouette calculation will use non-sieved frames ONLY.\n");
    }
  }
  if (!halffile_.empty()) {
    mprintf("\tSummary comparing parts of trajectory data for clusters will be written to %s\n",
            halffile_.c_str());
    if (!splitFrames_.empty()) {
      mprintf("\t\tFrames will be split at:");
      for (std::vector<int>::const_iterator f = splitFrames_.begin(); f != splitFrames_.end(); ++f)
        mprintf(" %i", *f);
      mprintf("\n");
    } else
      mprintf("\t\tFrames will be split at the halfway point.\n");
  }
  mprintf("\tRepresentative frames will be chosen by");
  switch (bestRep_) {
    case CUMULATIVE: mprintf(" lowest cumulative distance to all other frames.\n"); break;
    case CENTROID  : mprintf(" closest distance to cluster centroid.\n"); break;
    case CUMULATIVE_NOSIEVE:
      mprintf(" lowest cumulative distance to all other frames (ignore sieved frames).\n");
      break;
  }
  if (nRepsToSave_ > 1)
    mprintf("\tThe top %i representative frames will be determined.\n", nRepsToSave_);
  if (!clusterfile_.empty())
    mprintf("\tCluster trajectories will be written to %s, format %s\n",
            clusterfile_.c_str(), TrajectoryFile::FormatString(clusterfmt_));
  if (!singlerepfile_.empty())
    mprintf("\tCluster representatives will be written to 1 traj (%s), format %s\n",
            singlerepfile_.c_str(), TrajectoryFile::FormatString(singlerepfmt_));
  if (!reptrajfile_.empty()) {
    mprintf("\tCluster representatives will be written to separate trajectories,\n");
    mprintf("\t\tprefix (%s), format %s",reptrajfile_.c_str(), 
            TrajectoryFile::FormatString(reptrajfmt_));
    if (writeRepFrameNum_) mprintf(", with frame #s");
    mprintf("\n");
  }
  if (!avgfile_.empty())
    mprintf("\tAverage structures for clusters will be written to %s, format %s\n",
            avgfile_.c_str(), TrajectoryFile::FormatString(avgfmt_));
  if (!refs_.empty())
    mprintf("\tClusters will be identified with loaded reference structures if RMSD\n"
            "\t  (mask '%s') to representative frame is < %g Ang.\n", refmaskexpr_.c_str(),refCut_);
  if (drawGraph_ > 0)
    mprintf("\tEXPERIMENTAL: Force-directed graph will be drawn from pairwise distances.\n"
            "\t              Max iterations= %i, min tolerance= %g\n",
                             draw_maxit_, draw_tol_);

  return Analysis::OK;
}

/** This is where the clustering is actually performed. First the distances
  * between each frame are calculated. Then the clustering routine is called.
  */
// TODO: Need to update save to indicate distance type
// NOTE: Should distances be saved only if load_pair?
Analysis::RetType Analysis_Clustering::Analyze() {
  Timer cluster_setup;
  Timer cluster_pairwise;
  Timer cluster_cluster;
  Timer cluster_post;
  Timer cluster_total;
  cluster_total.Start();
  mprintf("\tStarting clustering.\n");
  cluster_setup.Start();
  // If no dataset specified, use COORDS
  if (cluster_dataset_.empty()) {
    if (coords_ == 0) {
      mprinterr("Error: No data to cluster on.\n");
      return Analysis::ERR;
    }
    cluster_dataset_.push_back( (DataSet*)coords_ );
  }
  // Test that cluster data set contains data
  unsigned int clusterDataSetSize = cluster_dataset_[0]->Size();
  if (clusterDataSetSize < 1) {
    mprinterr("Error: cluster data set %s does not contain data.\n", 
              cluster_dataset_[0]->legend());
    return Analysis::ERR;
  }
  // If more than one data set, make sure they are all the same size.
  for (ClusterDist::DsArray::iterator ds = cluster_dataset_.begin();
                                      ds != cluster_dataset_.end(); ++ds)
  {
    if ((*ds)->Size() != clusterDataSetSize) {
      mprinterr("Error: data set '%s' size (%zu) != first data set '%s' size (%u)\n",
                (*ds)->legend(), (*ds)->Size(), 
                cluster_dataset_[0]->legend(), clusterDataSetSize);
      return Analysis::ERR;
    }
  }
  // If no coordinates were specified, disable coordinate output types
  bool has_coords = true;
  if (coords_ == 0 || coords_->Size() < 1) {
    mprintf("Warning: No coordinates or associated coordinate data set is empty.\n"
            "Warning: Disabling coordinate output.\n");
    has_coords = false;
  }
  cluster_setup.Stop();

  // Set up cluster distance calculation
  if (CList_->SetupCdist( cluster_dataset_, metric_, nofitrms_, useMass_, maskexpr_ ))
    return Analysis::ERR;

  // Check or calculate pairwise distances between frames
  cluster_pairwise.Start();
  // Do some sanity checking first
  if (pw_dist_ == 0) {
    mprinterr("Internal Error: Empty cluster matrix.\n");
    return Analysis::ERR;
  }
  if (pw_dist_->Group() != DataSet::CLUSTERMATRIX) {
    mprinterr("Internal Error: Wrong cluster matrix type.\n");
    return Analysis::ERR;
  }
  if (pw_dist_->Size() > 0) {
    // Check if existing pw dist matrix matches expected size. If not, need to
    // allocate a new data set and recalculate.
    int sval = sieve_;
    if (sval < 0) sval = -sval;
    unsigned int expected_nrows = cluster_dataset_[0]->Size() / (unsigned int)sval;
    if ( (cluster_dataset_[0]->Size() % (unsigned int)sval) != 0 )
      expected_nrows++;
    if ( ((DataSet_Cmatrix*)pw_dist_)->Nrows() != expected_nrows ) {
      mprintf("Warning: ClusterMatrix has %zu rows, expected %zu; recalculating matrix.\n",
              ((DataSet_Cmatrix*)pw_dist_)->Nrows(), expected_nrows);
      pw_dist_ = masterDSL_->AddSet(DataSet::CMATRIX, "", "CMATRIX");
      if (pw_dist_ == 0) return Analysis::ERR;
    }
  }
  if (CList_->CalcFrameDistances( pw_dist_, cluster_dataset_, sieve_, sieveSeed_ ))
    return Analysis::ERR;
  // If we want to save the pairwise distances add to file now
  if (pwd_file_ != 0)
    pwd_file_->AddDataSet( pw_dist_ );
  cluster_pairwise.Stop();

  // Cluster
  cluster_cluster.Start();
  CList_->Cluster();
  cluster_cluster.Stop();
  cluster_post.Start();
  Timer cluster_post_renumber;
  Timer cluster_post_bestrep;
  Timer cluster_post_info;
  Timer cluster_post_summary;
  Timer cluster_post_coords;
  if (CList_->Nclusters() > 0) {
    // Sort clusters and renumber; also finds centroids. If sieving,
    // add remaining frames.
    cluster_post_renumber.Start();
    CList_->Renumber( (sieve_ != 1) );
    cluster_post_renumber.Stop();
    // Find best representative frames for each cluster.
    cluster_post_bestrep.Start();
    switch (bestRep_) {
      case CUMULATIVE: CList_->FindBestRepFrames_CumulativeDist(nRepsToSave_); break;
      case CENTROID  : CList_->FindBestRepFrames_Centroid(nRepsToSave_); break;
      case CUMULATIVE_NOSIEVE: CList_->FindBestRepFrames_NoSieve_CumulativeDist(nRepsToSave_); break;
    }
    cluster_post_bestrep.Stop();
    // DEBUG
    if (debug_ > 0) {
      mprintf("\nFINAL CLUSTERS:\n");
      CList_->PrintClusters();
    }
    // Attempt to assign reference names to clusters if any specified.
    if (!refs_.empty()) {
      if (has_coords)
        AssignRefsToClusters( *CList_ );
      else
        mprintf("Warning: References were specified but no COORDS. Cannot assign ref names.\n");
    }

    // Print ptraj-like cluster info.
    // If no filename is written and no noinfo, some info will still be written to STDOUT
    if (!suppressInfo_) {
      cluster_post_info.Start();
      CList_->PrintClustersToFile(clusterinfo_);
      cluster_post_info.Stop();
    }

    // Calculate cluster silhouette
    if (!sil_file_.empty())
      CList_->CalcSilhouette( sil_file_, includeSieveInCalc_ );

    // Print a summary of clusters
    if (!summaryfile_.empty()) {
      cluster_post_summary.Start();
      CList_->Summary(summaryfile_, includeSieveInCalc_);
      cluster_post_summary.Stop();
    }

    // Print a summary comparing first half to second half of data for clusters
    if (!halffile_.empty()) {
      // If no split frames were specified, use halfway point.
      if (splitFrames_.empty())
        splitFrames_.push_back( clusterDataSetSize / 2 );
      // Check that none of the split values are invalid.
      std::vector<int> actualSplitFrames;
      for (std::vector<int>::const_iterator f = splitFrames_.begin();
                                            f != splitFrames_.end(); ++f)
        if ( *f < 1 || *f >= (int)clusterDataSetSize )
          mprintf("Warning: split frame %i is out of bounds; ignoring.\n", *f);
        else
          actualSplitFrames.push_back( *f );
      CList_->Summary_Part(halffile_, actualSplitFrames);
    }

    // Create cluster v time data from clusters.
    CreateCnumvtime( *CList_, clusterDataSetSize );

    // TEST: Draw graph based on point distances
    if (drawGraph_ > 0)
     CList_->DrawGraph( drawGraph_ == 2, cnumvtime_, draw_tol_, draw_maxit_ );

    // Create # clusters seen v time data.
    if (clustersVtime_ != 0)
      NclustersObserved( *CList_, clusterDataSetSize );

    // Create cluster pop v time plots
    if (cpopvtimefile_ != 0)
      CreateCpopvtime( *CList_, clusterDataSetSize );

    // Create cluster lifetime DataSets
    if (calc_lifetimes_)
      ClusterLifetimes( *CList_, clusterDataSetSize );

    // Change cluster num v time to grace-compatible colors if specified.
    if (grace_color_) {
      DataSet_integer& cnum_temp = static_cast<DataSet_integer&>( *cnumvtime_ );
      for (DataSet_integer::iterator ival = cnum_temp.begin();
                                     ival != cnum_temp.end(); ++ival)
      {
        *ival += 1;
        if (*ival > 15) *ival = 15;
      }
    }
    // Coordinate output.
    if (has_coords) {
      cluster_post_coords.Start();
      // Write clusters to trajectories
      if (!clusterfile_.empty())
        WriteClusterTraj( *CList_ ); 
      // Write all representative frames to a single traj
      if (!singlerepfile_.empty())
        WriteSingleRepTraj( *CList_ );
      // Write all representative frames to separate trajs
      if (!reptrajfile_.empty())
        WriteRepTraj( *CList_ );
      // Write average structures for each cluster to separate files.
      if (!avgfile_.empty())
        WriteAvgStruct( *CList_ );
      cluster_post_coords.Stop();
    }
  } else
    mprintf("\tNo clusters found.\n");
  cluster_post.Stop();
  cluster_total.Stop();
  // Timing data
  mprintf("\tCluster timing data:\n");
  cluster_setup.WriteTiming(1,    "  Cluster Init. :", cluster_total.Total());
  cluster_pairwise.WriteTiming(1, "  Pairwise Calc.:", cluster_total.Total());
  cluster_cluster.WriteTiming(1,  "  Clustering    :", cluster_total.Total());
# ifdef TIMER
  CList_->Timing( cluster_cluster.Total() );
# endif
  cluster_post.WriteTiming(1,     "  Cluster Post. :", cluster_total.Total());
  cluster_post_renumber.WriteTiming(2, "Cluster renumbering/sieve restore", cluster_post.Total());
  cluster_post_bestrep.WriteTiming(2, "Find best rep.", cluster_post.Total());
  cluster_post_info.WriteTiming(2, "Info calc", cluster_post.Total());
  cluster_post_summary.WriteTiming(2, "Summary calc", cluster_post.Total());
  cluster_post_coords.WriteTiming(2, "Coordinate writes", cluster_post.Total());
  cluster_total.WriteTiming(1,    "Total:");
  return Analysis::OK;
}

// -----------------------------------------------------------------------------
void Analysis_Clustering::AssignRefsToClusters( ClusterList& CList ) const {
  // Pre-center all reference coords at the origin. No need to store trans vectors.
  std::vector<Frame> refFrames;
  refFrames.reserve( refs_.size() );
  for (unsigned int idx = 0; idx != refs_.size(); idx++) {
    AtomMask rMask( refmaskexpr_ );
    DataSet_Coords_REF* REF_ds = (DataSet_Coords_REF*)refs_[idx];
    if ( REF_ds->Top().SetupIntegerMask( rMask, REF_ds->RefFrame() ) ) {
      mprintf("Warning: Could not set up mask for reference '%s'\n", REF_ds->legend());
      continue;
    }
    refFrames.push_back( Frame(REF_ds->RefFrame(), rMask) );
    refFrames.back().CenterOnOrigin( useMass_ );
  }
  // For each cluster, assign the reference name with the lowest RMSD
  // to the representative frame that is below the cutoff.
  AtomMask tMask( refmaskexpr_ );
  if (coords_->Top().SetupIntegerMask( tMask )) {
    mprinterr("Error: Could not set up mask for assigning references.\n");
    return;
  }
  Frame TGT( coords_->AllocateFrame(), tMask );
  unsigned int cidx = 0;
  for (ClusterList::cluster_it cluster = CList.begin();
                               cluster != CList.end(); ++cluster, ++cidx)
  {
    coords_->GetFrame( cluster->BestRepFrame(), TGT, tMask );
    double minRms = TGT.RMSD_CenteredRef( refFrames[0], useMass_ );
    unsigned int minIdx = 0;
    for (unsigned int idx = 1; idx < refs_.size(); idx++) {
      double rms = TGT.RMSD_CenteredRef( refFrames[idx], useMass_ );
      if (rms < minRms) {
        minRms = rms;
        minIdx = idx;
      }
    }
    if (minRms < refCut_) {
      //mprintf("DEBUG: Assigned cluster %i to reference \"%s\" (%g)\n", cidx,
      //        refs_[minIdx]->Meta().Name().c_str(), minRms);
      cluster->SetNameAndRms( refs_[minIdx]->Meta().Name(), minRms );
    } else {
      //mprintf("DEBUG: Cluster %i was closest to reference \"(%s)\" (%g)\n", cidx,
      //        refs_[minIdx]->Meta().Name().c_str(), minRms);
      cluster->SetNameAndRms( "(" + refs_[minIdx]->Meta().Name() + ")", minRms );
    }
  }
}

// -----------------------------------------------------------------------------
// Analysis_Clustering::CreateCnumvtime()
/** Put cluster number vs frame into dataset.  */
void Analysis_Clustering::CreateCnumvtime( ClusterList const& CList, unsigned int maxFrames ) {
  // FIXME:
  // Cast generic DataSet for cnumvtime back to integer dataset to 
  // access specific integer dataset functions for resizing and []
  // operator. Should this eventually be generic to all atomic DataSets? 
  DataSet_integer& cnum_temp = static_cast<DataSet_integer&>( *cnumvtime_ );
  cnum_temp.Resize( maxFrames );
  // Make all clusters start at -1. This way cluster algorithms that
  // have noise points (i.e. no cluster assigned) will be distinguished.
  std::fill(cnum_temp.begin(), cnum_temp.end(), -1);

  for (ClusterList::cluster_iterator C = CList.begincluster();
                                     C != CList.endcluster(); C++)
  {
    //mprinterr("Cluster %i:\n",CList->CurrentNum());
    int cnum = (*C).Num();
    // Loop over all frames in the cluster
    for (ClusterNode::frame_iterator frame = (*C).beginframe();
                                     frame != (*C).endframe(); frame++)
    {
      //mprinterr("%i,",*frame);
      cnum_temp[ *frame ] = cnum;
    }
    //mprinterr("\n");
    //break;
  }
}

// Analysis_Clustering::CreateCpopvtime()
// NOTE: Should not be called if cpopvtimefile is NULL
void Analysis_Clustering::CreateCpopvtime( ClusterList const& CList, unsigned int maxFrames ) {
  mprintf("\tCalculating cluster population vs time for each cluster.\n");
  // Set up output data sets
  std::vector<DataSet_float*> Cpop;
  MetaData md(cnumvtime_->Meta().Name(), "Pop");
  DataSet::SizeArray dsize(1, maxFrames);
  for (int cnum = 0; cnum < CList.Nclusters(); ++cnum) {
    md.SetIdx( cnum );
    DataSet_float* ds = (DataSet_float*)masterDSL_->AddSet( DataSet::FLOAT, md );
    if (ds == 0) {
      mprinterr("Error: Could not allocate cluster pop v time DataSet.\n");
      return;
    }
    ds->Allocate( dsize );
    Cpop.push_back( ds );
    // SANITY CHECK
    if (cpopvtimefile_ != 0)
      cpopvtimefile_->AddDataSet( ds );
  }
  int cnum = 0;
  // Loop over all clusters
  for (ClusterList::cluster_iterator C = CList.begincluster(); C != CList.endcluster(); ++C, ++cnum)
  {
    DataSet_float& pvt = static_cast<DataSet_float&>( *(Cpop[cnum]) );
    float pop = 0.0;
    // Loop over all frames in cluster
    for (ClusterNode::frame_iterator f = C->beginframe(); f != C->endframe(); ++f)
    {
      if (*f > (int)pvt.Size())
        pvt.Resize( *f, pop );
      pop = pop + 1.0;
      pvt[*f] = pop;
    }
    // Ensure pop v time set is maxFrames long
    if (pvt.Size() < maxFrames)
      pvt.Resize( maxFrames, pop );
    // Normalization
    if (norm_pop_ == CLUSTERPOP) {
      float norm = 1.0 / (float)C->Nframes();
      for (unsigned int frm = 0; frm < maxFrames; ++frm)
        pvt[frm] = pvt[frm] * norm;
    } else if (norm_pop_ == FRAME) {
      float norm = 1.0;
      for (unsigned int frm = 0; frm < maxFrames; ++frm)
      {
        pvt[frm] = pvt[frm] / norm;
        norm = norm + 1.0;
      }
    }
  }
}

// Analysis_Clustering::ClusterLifetimes()
void Analysis_Clustering::ClusterLifetimes( ClusterList const& CList, unsigned int maxFrames ) {
  // Set up output data sets. TODO: use ChildDSL
  std::vector<DataSet_integer*> DSL;
  MetaData md(cnumvtime_->Meta().Name(), "Lifetime");
  for (int cnum = 0; cnum < CList.Nclusters(); ++cnum) {
    md.SetIdx( cnum ); 
    DSL.push_back((DataSet_integer*) masterDSL_->AddSet( DataSet::INTEGER, md));
    if (DSL.back() == 0) {
      mprinterr("Error: Could not allocate cluster lifetime DataSet.\n");
      return;
    }
    DSL.back()->Resize( maxFrames );
  }
  // For each frame, assign cluster frame belongs to 1.
  DataSet_integer const& cnum_temp = static_cast<DataSet_integer const&>(*cnumvtime_);
  for (unsigned int frame = 0; frame < maxFrames; ++frame) {
    int cluster_num = cnum_temp[frame];
    // Noise points are -1
    if (cluster_num > -1)
      (*DSL[ cluster_num ])[ frame ] = 1;
  }
}

/** Determine how many different clusters are observed within a given time
  * window.
  */
void Analysis_Clustering::NclustersObserved( ClusterList const& CList, unsigned int maxFrames ) {
  DataSet_integer const& CVT = static_cast<DataSet_integer const&>( *cnumvtime_ );
  if (CVT.Size() < 1 || CList.Nclusters() < 1) return;
  int dataIdx = 0;
  // True if cluster was observed during window
  std::vector<bool> observed( CList.Nclusters(), false );
  for (unsigned int frame = 0; frame < maxFrames; frame++) {
    if (CVT[frame] != -1)
      observed[ CVT[frame] ] = true;
    if ( ((frame+1) % windowSize_) == 0 ) {
      // Count # observed clusters
      int nClustersObserved = 0;
      for (std::vector<bool>::iterator ob = observed.begin(); ob != observed.end(); ++ob)
        if ( *ob ) {
          ++nClustersObserved;
          *ob = false;
        }
      //mprintf("DEBUG: WINDOW at frame %i; %i clusters observed\n", frame+1, nClustersObserved);
      clustersVtime_->Add( dataIdx++, &nClustersObserved );
    }
  }
/*
  int currentCluster = CVT[0];
  int nClustersObserved = 1;
  for (int frame = 1; frame < maxFrames; frame++) {
    // Do not count noise as a cluster.
    if (CVT[frame] != currentCluster && CVT[frame] != -1) {
      ++nClustersObserved;
      currentCluster = CVT[frame];
    }
    mprintf("DEBUG: %i %i\n", frame+1, nClustersObserved);
    if ( ((frame+1) % windowSize_) == 0 ) {
      mprintf("DEBUG: WINDOW\n");
      clustersVtime_->Add( dataIdx++, &nClustersObserved );
      nClustersObserved = 1;
    }
  }
*/
  clustersVtime_->SetDim(Dimension::X, Dimension(windowSize_, windowSize_, "Frame"));
} 

// ---------- Cluster Coordinate Output Routines -------------------------------
// Analysis_Clustering::WriteClusterTraj()
/** Write frames in each cluster to a trajectory file.  */
void Analysis_Clustering::WriteClusterTraj( ClusterList const& CList ) {
  Topology* clusterparm = coords_->TopPtr();
  // Loop over all clusters
  for (ClusterList::cluster_iterator C = CList.begincluster();
                                     C != CList.endcluster(); ++C)
  {
    // Create filename based on cluster number.
    int cnum = C->Num();
    std::string cfilename =  clusterfile_ + ".c" + integerToString( cnum );
    // Set up trajectory file 
    Trajout_Single clusterout;
    if (clusterout.PrepareTrajWrite(cfilename, ArgList(), clusterparm,
                                    coords_->CoordsInfo(), C->Nframes(),
                                    clusterfmt_)) 
    {
      mprinterr("Error: Could not set up cluster trajectory %s for write.\n",
                cfilename.c_str());
      return;
    }
    // Loop over all frames in cluster
    int set = 0;
    Frame clusterframe = coords_->AllocateFrame();
    for (ClusterNode::frame_iterator fnum = C->beginframe();
                                     fnum != C->endframe(); ++fnum)
    {
      coords_->GetFrame( *fnum, clusterframe );
      clusterout.WriteSingle(set++, clusterframe);
    }
    // Close traj
    clusterout.EndTraj();
  }
}

// Analysis_Clustering::WriteAvgStruct()
void Analysis_Clustering::WriteAvgStruct( ClusterList const& CList ) {
  Topology avgparm = coords_->Top();
  // Get extension for representative frame format 
  std::string tmpExt = TrajectoryFile::WriteFormatExtension(avgfmt_);
  // Loop over all clusters
  for (ClusterList::cluster_iterator C = CList.begincluster();
                                     C != CList.endcluster(); ++C)
  {
    // Create filename based on cluster number.
    int cnum = C->Num();
    std::string cfilename = avgfile_ + ".c" + integerToString( cnum ) + tmpExt;
    // Set up trajectory file
    Trajout_Single clusterout; // FIXME CoordinateInfo OK for just coords?
    if (clusterout.PrepareTrajWrite(cfilename, ArgList(), &avgparm,
                                    CoordinateInfo(), 1, avgfmt_))
    {
      mprinterr("Error: Could not set up cluster average file %s for write.\n",
                cfilename.c_str());
      return;
    }
    // Get rep frame for rms fitting.
    Frame repframe = coords_->AllocateFrame();
    coords_->GetFrame( C->BestRepFrame(), repframe );
    Vec3 reftrans = repframe.CenterOnOrigin(false);
    // Loop over all frames in cluster
    Frame clusterframe = coords_->AllocateFrame();
    Frame avgframe = clusterframe;
    avgframe.ZeroCoords();
    for (ClusterNode::frame_iterator fnum = C->beginframe();
                                     fnum != C->endframe(); ++fnum)
    {
      coords_->GetFrame( *fnum, clusterframe );
      clusterframe.RMSD_FitToRef( repframe, reftrans );
      avgframe += clusterframe;
    }
    avgframe.Divide( (double)C->Nframes() );
    clusterout.WriteSingle(0, avgframe);
    clusterout.EndTraj();
  }
}
 
// Analysis_Clustering::WriteSingleRepTraj()
/** Write representative frame of each cluster to a trajectory file.  */
void Analysis_Clustering::WriteSingleRepTraj( ClusterList const& CList ) {
  Trajout_Single clusterout;
  // Set up trajectory file. Use parm from COORDS DataSet. 
  Topology *clusterparm = coords_->TopPtr();
  if (clusterout.PrepareTrajWrite(singlerepfile_, ArgList(), clusterparm,
                                  coords_->CoordsInfo(), CList.Nclusters() * nRepsToSave_,
                                  singlerepfmt_)) 
  {
    mprinterr("Error: Could not set up single trajectory for represenatatives %s for write.\n",
                singlerepfile_.c_str());
     return;
  }
  // Set up frame to hold cluster rep coords. 
  Frame clusterframe = coords_->AllocateFrame();
  int framecounter = 0;
  // Write rep frames from all clusters.
  for (ClusterList::cluster_iterator cluster = CList.begincluster(); 
                                     cluster != CList.endcluster(); ++cluster) 
  {
    for (ClusterNode::RepPairArray::const_iterator rep = cluster->BestReps().begin();
                                                   rep != cluster->BestReps().end(); ++rep)
    {
      coords_->GetFrame( rep->first, clusterframe );
      clusterout.WriteSingle(framecounter++, clusterframe);
    }
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
  std::string tmpExt = TrajectoryFile::WriteFormatExtension(reptrajfmt_);
  // Use Topology from COORDS DataSet to set up input frame
  Topology* clusterparm = coords_->TopPtr();
  Frame clusterframe = coords_->AllocateFrame();
  // Loop over all clusters
  for (ClusterList::cluster_iterator C = CList.begincluster();
                                     C != CList.endcluster(); ++C)
  {
    if (writeRepFrameNum_) {
      // Each rep from cluster to separate file.
      for (ClusterNode::RepPairArray::const_iterator rep = C->BestReps().begin();
                                                     rep != C->BestReps().end(); ++rep)
      {
        Trajout_Single clusterout;
        // Get best rep frame # 
        int framenum = rep->first;
        // Create filename based on cluster number and frame #
        std::string cfilename = reptrajfile_ + ".c" + integerToString(C->Num()) +
                                ("." + integerToString(framenum+1)) + tmpExt;
        // Set up trajectory file.
        if (clusterout.PrepareTrajWrite(cfilename, ArgList(), clusterparm,
                                        coords_->CoordsInfo(), 1, reptrajfmt_))
        {
          mprinterr("Error: Could not set up representative trajectory file %s for write.\n",
                    cfilename.c_str());
          return;
        }
        // Write cluster rep frame
        coords_->GetFrame( framenum, clusterframe );
        clusterout.WriteSingle(framenum, clusterframe);
        // Close traj
        clusterout.EndTraj();
      }
    } else {
      // Each rep from cluster to single file.
      Trajout_Single clusterout;
      // Create filename based on cluster number
      std::string cfilename = reptrajfile_ + ".c" + integerToString(C->Num()) + tmpExt;
      // Set up trajectory file.
      if (clusterout.PrepareTrajWrite(cfilename, ArgList(), clusterparm,
                                      coords_->CoordsInfo(), nRepsToSave_, reptrajfmt_))
      {
        mprinterr("Error: Could not set up representative trajectory file %s for write.\n",
                  cfilename.c_str());
        return;
      }
      int framecounter = 0;
      for (ClusterNode::RepPairArray::const_iterator rep = C->BestReps().begin();
                                                     rep != C->BestReps().end(); ++rep)
      {
        // Write cluster rep frame
        coords_->GetFrame( rep->first, clusterframe );
        clusterout.WriteSingle( framecounter++, clusterframe );
      }
      // Close traj
      clusterout.EndTraj();
    }
  }
}

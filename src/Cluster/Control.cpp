#include "Control.h"
#include "../CpptrajStdio.h"
#include "../DataSet_Coords.h"
#include "Output.h"
#include "../DataFile.h" // For loading pairwise cache
// Metric classes
#include "Metric_RMS.h"
// Algorithms
#include "Algorithm_HierAgglo.h"
#include "Algorithm_DBscan.h"

Cpptraj::Cluster::Control::Control() :
  metric_(0),
  algorithm_(0),
  verbose_(0),
  frameSelect_(UNSPECIFIED),
  sieve_(1),
  sieveSeed_(-1),
  sieveRestore_(NO_RESTORE),
  restoreEpsilon_(0.0),
  includeSieveInCalc_(false),
  bestRep_(BestReps::NO_REPS),
  nRepsToSave_(1),
  suppressInfo_(false),
  cnumvtime_(0),
  cpopvtimefile_(0),
  norm_pop_(Node::NONE)
{}

Cpptraj::Cluster::Control::~Control() {
  if (algorithm_ != 0) delete algorithm_;
  if (metric_ != 0   ) delete metric_;
}

// -----------------------------------------------------------------------------
/** The default pairwise cache file name. */
const char* Cpptraj::Cluster::Control::DEFAULT_PAIRDIST_NAME_ = "CpptrajPairDist";

/** The default pairwise distance file type. */
DataFile::DataFormatType Cpptraj::Cluster::Control::DEFAULT_PAIRDIST_TYPE_ =
# ifdef BINTRAJ
  DataFile::CMATRIX_NETCDF;
# else
  DataFile::CMATRIX_BINARY;
# endif

const char* Cpptraj::Cluster::Control::PairwiseArgs_ =
  "[pairdist <name>] [loadpairdist] [savepairdist] [pairwisecache {mem|disk|none}]";

/** Set up PairwiseMatrix from arguments. */
int Cpptraj::Cluster::Control::AllocatePairwise(ArgList& analyzeArgs, DataSetList& DSL,
                                                DataFileList& DFL)
{
  if (metric_ == 0) {
    mprinterr("Internal Error: AllocatePairwise(): Metric is null.\n");
    return 1;
  }

  // Determine if we are saving/loading pairwise distances
  std::string pairdistname = analyzeArgs.GetStringKey("pairdist");
  bool load_pair = analyzeArgs.hasKey("loadpairdist");
  bool save_pair = analyzeArgs.hasKey("savepairdist");
  // Check if we need to set a default file name
  std::string fname;
  if (pairdistname.empty())
    fname = DEFAULT_PAIRDIST_NAME_;
  else {
    fname = pairdistname;
    // To remain backwards compatible, assume we want to load if
    // a pairdist name was specified.
    if (!load_pair && !save_pair) {
      mprintf("Warning: 'pairdist' specified but 'loadpairdist'/'savepairdist' not specified."
              "Warning: Assuming 'loadpairdist'.\n");
      load_pair = true;
    }
  }

  cache_ = 0;
  if (load_pair) {
    DataSetList selected;
    // If specified, check if the pairdistname refers to a DataSet
    if (!pairdistname.empty())
      selected = DSL.SelectGroupSets( pairdistname, DataSet::PWCACHE );
    if (selected.empty()) {
      // No DataSet; check if we can load from a file.
      if (File::Exists( fname )) {
        DataFile dfIn;
        if (dfIn.ReadDataIn( fname, ArgList(), DSL )) return 1;
        DataSet* ds = DSL.GetDataSet( fname );
        if (ds == 0) return 1;
        if (ds->Group() != DataSet::PWCACHE) {
          mprinterr("Internal Error: AllocatePairwise(): Set is not a pairwise cache.\n");
          return 1;
        }
        cache_ = (DataSet_PairwiseCache*)ds;
        mprintf("DEBUG: Loaded cache set '%s' from file: %s\n",
                cache_->legend(), dfIn.DataFilename().full());
        // No need to save pw cache if we just loaded one.
        if (save_pair) {
          mprintf("Warning: 'savepairdist' specified but pairwise cache loaded from file.\n"
                  "Warning: Disabling 'savepairdist'.\n");
          save_pair = false;
        }
      }
    } else {
      if (selected.size() > 1)
        mprintf("Warning: '%s' matches multiple sets; only using '%s'\n",
                pairdistname.c_str(), selected[0]->legend());
      cache_ = (DataSet_PairwiseCache*)selected[0];
    }
    if (cache_ == 0)
      mprintf("Warning: 'loadpairdist' specified but cache set/file not found.\n");
  } // END if load_pair

  if (cache_ == 0) {
    // Process DataSet type arguments
    DataSet::DataType pw_type = DataSet::PMATRIX_MEM;
    std::string pw_typeString = analyzeArgs.GetStringKey("pairwisecache");
    if (!pw_typeString.empty()) {
      if (pw_typeString == "mem")
        pw_type = DataSet::PMATRIX_MEM; 
      else if (pw_typeString == "disk")
        pw_type = DataSet::PMATRIX_NC;
      else if (pw_typeString == "none")
        pw_type = DataSet::UNKNOWN_DATA;
      else {
        mprinterr("Error: Unrecognized option for 'pairwisecache' ('%s')\n", pw_typeString.c_str());
        return 1;
      }
    }
    // Allocate cache if necessary
    if (pw_type != DataSet::UNKNOWN_DATA) {
      MetaData meta( pairdistname );
      // Cache-specific setup.
      if (pw_type == DataSet::PMATRIX_NC)
        meta.SetFileName( fname );
      cache_ = (DataSet_PairwiseCache*)DSL.AddSet( pw_type, meta, "CMATRIX" );
      if (cache_ == 0) {
        mprinterr("Error: Could not allocate pairwise cache.\n");
        return 1;
      }
      mprintf("DEBUG: Allocated pairwise distance cache: %s\n", cache_->legend());
    }
  }

  // Setup pairwise matrix
  if (pmatrix_.Setup(metric_, cache_)) return 1;

  if (save_pair) {
    if (cache_ == 0) {
      mprintf("Warning: Not caching distances; ignoring 'savepairdist'\n");
    } else {
      DataFile::DataFormatType pw_file_type = DataFile::UNKNOWN_DATA;
      if (pairdistname.empty())
        pw_file_type = DEFAULT_PAIRDIST_TYPE_;
      // TODO enable saving for other set types?
      if (cache_->Type() == DataSet::PMATRIX_MEM) {
        DataFile* pwd_file = DFL.AddDataFile( fname, pw_file_type, ArgList() );
        if (pwd_file == 0) return 1;
        pwd_file->AddDataSet( cache_ );
        mprintf("DEBUG: Saving pw distance cache '%s' to file '%s'\n", cache_->legend(),
                pwd_file->DataFilename().full());
      } else
        mprintf("Warning: Can only save pairwise distance cache for in-memory caches.\n");
    }
  }

  return 0;
}

// -----------------------------------------------------------------------------
/** \return Pointer to Algorithm of given type. */
Cpptraj::Cluster::Algorithm* Cpptraj::Cluster::Control::AllocateAlgorithm(Algorithm::AType atype)
{
  Algorithm* alg = 0;
  switch (atype) {
    case Algorithm::HIERAGGLO : alg = new Algorithm_HierAgglo(); break;
    case Algorithm::DBSCAN    : alg = new Algorithm_DBscan(); break;
    default : mprinterr("Error: Unhandled Algorithm in AllocateAlgorithm.\n");
  }
  return alg;
}

const char* Cpptraj::Cluster::Control::AlgorithmArgs_ =
  "{hieragglo|dbscan|kmeans|dpeaks}";

/** Set up Algorithm from keyword + arguments. */
int Cpptraj::Cluster::Control::AllocateAlgorithm(ArgList& analyzeArgs) {
  Algorithm::AType atype;
  if      (analyzeArgs.hasKey("hieragglo")) atype = Algorithm::HIERAGGLO;
  else if (analyzeArgs.hasKey("dbscan"   )) atype = Algorithm::DBSCAN;
  else if (analyzeArgs.hasKey("kmeans") ||
           analyzeArgs.hasKey("means")    ) atype = Algorithm::KMEANS;
  else if (analyzeArgs.hasKey("dpeaks"   )) atype = Algorithm::DPEAKS;
  else {
    mprintf("Warning: No clustering algorithm specified; defaulting to 'hieragglo'\n");
    atype = Algorithm::HIERAGGLO;
  }
  if (algorithm_ != 0) delete algorithm_;
  algorithm_ = AllocateAlgorithm( atype );
  if (algorithm_ == 0) return 1;
  if (algorithm_->Setup( analyzeArgs )) return 1;
  return 0;
}

// -----------------------------------------------------------------------------
/** /return Pointer to Metric of given type. */
Cpptraj::Cluster::Metric* Cpptraj::Cluster::Control::AllocateMetric(Metric::Type mtype)
{
  Metric* met = 0;
  switch (mtype) {
    case Metric::RMS : met = new Metric_RMS(); break;
    default: mprinterr("Error: Unhandled Metric in AllocateMetric.\n");
  }
  return met;
}

const char* Cpptraj::Cluster::Control::CoordsDataSetArgs_ =
  "{dme|rms|srmsd} [mass] [nofit] [<mask>]";

/** Set up clustering for a COORDS DataSet.*/
int Cpptraj::Cluster::Control::SetupForCoordsDataSet(DataSet_Coords* ds,
                                                     ArgList& analyzeArgs,
                                                     DataSetList& DSL,
                                                     DataFileList& DFL,
                                                     int verboseIn)
{
  verbose_ = verboseIn;
  if (ds == 0) {
    mprinterr("Internal Error: Control::SetupForCoordsDataSet() called with null DataSet.\n");
    return 1;
  }
  // Determine Metric. Valid ones for COORDS are RMS, DME, SRMSD
  int usedme = (int)analyzeArgs.hasKey("dme");
  int userms = (int)analyzeArgs.hasKey("rms");
  int usesrms = (int)analyzeArgs.hasKey("srmsd");
  if (usedme + userms + usesrms > 1) {
    mprinterr("Error: Specify either 'dme', 'rms', or 'srmsd'.\n");
    return 1;
  }
  Metric::Type mtype = Metric::RMS; // Default
  if      (usedme)  mtype = Metric::DME;
  else if (userms)  mtype = Metric::RMS;
  else if (usesrms) mtype = Metric::SRMSD;
  if (metric_ != 0) delete metric_;
  metric_ = 0;
  metric_ = AllocateMetric( mtype );
  if (metric_ == 0) return 1;
  // Metric setup.
  bool useMass = analyzeArgs.hasKey("mass");
  bool nofit   = analyzeArgs.hasKey("nofit");
  // Get the mask string 
  std::string maskExpr = analyzeArgs.GetMaskNext();
  /*if (!refs_.empty() && refmaskexpr_.empty()) { TODO enable
    refmaskexpr_ = maskexpr_;
    if (refmaskexpr_.empty()) {
      refmaskexpr_.assign("!@H=");
      mprintf("Warning: 'assignrefs' specified but no 'refmask' given.\n"
              "Warning:   Using default mask expression: '%s'\n", refmaskexpr_.c_str());
    }
  }*/

  int err = 0;
  switch (mtype) {
    case Metric::RMS :
      err = ((Metric_RMS*)metric_)->Init(ds, AtomMask(maskExpr), nofit, useMass); break;
    default:
      mprinterr("Error: Unhandled Metric setup.\n");
      err = 1;
  }
  if (err != 0) {
    mprinterr("Error: Metric setup failed.\n");
    return 1;
  }

  return Common(analyzeArgs, DSL, DFL);
}

// -----------------------------------------------------------------------------
const char* Cpptraj::Cluster::Control::CommonArgs_ =
  "[sieve <#> [sieveseed <#>] [random] [includesieveincalc] [sievetoframe]] "
  "[bestrep {cumulative|centroid|cumulative_nosieve} [savenreps <#>]] "
  "[noinfo|info <file>] [summary <file>] [sil <prefix>] "
  "[cpopvtime <file> [{normpop|normframe}]] "
  "[out <cnumvtime file>] [<ds name>] ";

/** Common setup. */
int Cpptraj::Cluster::Control::Common(ArgList& analyzeArgs, DataSetList& DSL, DataFileList& DFL)
{
  clusters_.SetDebug( verbose_ );

  // Allocate PairwiseMatrix (and optionally a cache). Metric must already be set up.
  if (AllocatePairwise( analyzeArgs, DSL, DFL )) {
    mprinterr("Error: PairwiseMatrix setup failed.\n");
    return 1;
  }

  // Allocate algorithm
  if (AllocateAlgorithm( analyzeArgs )) {
    mprinterr("Error: Algorithm setup failed.\n");
    return 1;
  }
  algorithm_->SetDebug( verbose_ );

  // Sieve options
  sieveSeed_ = analyzeArgs.getKeyInt("sieveseed", -1);
  sieve_ = analyzeArgs.getKeyInt("sieve", 1);
  if (sieve_ < 1) {
    mprinterr("Error: 'sieve <#>' must be >= 1 (%i)\n", sieve_);
    return 1;
  }
  if (analyzeArgs.hasKey("random") && sieve_ > 1)
    sieve_ = -sieve_; // negative # indicates random sieve

  // TODO incorporate with cumulative_nosieve? Or keep granular?
  includeSieveInCalc_ = analyzeArgs.hasKey("includesieveincalc");
  if (includeSieveInCalc_)
    mprintf("Warning: 'includesieveincalc' may be very slow.\n");

  // Determine how frames to cluster will be chosen
  if (frameSelect_ == UNSPECIFIED) {
    // If no other frame selection option like sieve provided and an already
    // set up cache is present, use the cached frames.
    if (sieve_ == 1 && pmatrix_.HasCache() && pmatrix_.Cache().Size() > 0)
    {
      frameSelect_ = FROM_CACHE;
      if (pmatrix_.Cache().SieveVal() != 1) {
        sieve_ = pmatrix_.Cache().SieveVal();
        mprintf("Warning: No sieve specified; using sieve value from cache '%s': %i\n",
                pmatrix_.Cache().legend(), sieve_);
      }
    }
  }

  // Choose sieve restore option based on algorithm and frames to cluster.
  if (sieve_ != 1)
  {
    sieveRestore_ = CLOSEST_CENTROID;
    if (algorithm_->Type() == Algorithm::DBSCAN) {
      if (!analyzeArgs.hasKey("sievetoframe"))
        sieveRestore_ = EPSILON_CENTROID;
      else
        sieveRestore_ = EPSILON_FRAME;
      restoreEpsilon_ = ((Algorithm_DBscan*)algorithm_)->Epsilon();
    }
  }

  // Best rep options
  std::string bestRepStr = analyzeArgs.GetStringKey("bestrep");
  if (bestRepStr.empty()) {
    // For sieving, cumulative can get very expensive. Default to centroid.
    if (sieve_ != 1)
      bestRep_ = BestReps::CENTROID;
    else
      bestRep_ = BestReps::CUMULATIVE;
  } else {
    if (bestRepStr == "cumulative")
      bestRep_ = BestReps::CUMULATIVE;
    else if (bestRepStr == "centroid")
      bestRep_ = BestReps::CENTROID;
    else if (bestRepStr == "cumulative_nosieve")
      bestRep_ = BestReps::CUMULATIVE_NOSIEVE;
    else {
      mprinterr("Error: Invalid 'bestRep' option (%s)\n", bestRepStr.c_str());
      return 1;
    }
  }
  nRepsToSave_ = analyzeArgs.getKeyInt("savenreps", 1);
  if (nRepsToSave_ < 1) {
    mprinterr("Error: 'savenreps' must be > 0\n");
    return 1;
  }

  // Cluster info output
  suppressInfo_ = analyzeArgs.hasKey("noinfo");
  if (!suppressInfo_)
    clusterinfo_ = analyzeArgs.GetStringKey("info");

  // Cluster summary output
  summaryfile_ = analyzeArgs.GetStringKey("summary");

  // Cluster silhouette output
  sil_file_ = analyzeArgs.GetStringKey("sil");

  // Cluster pop v time output
  cpopvtimefile_ = DFL.AddDataFile(analyzeArgs.GetStringKey("cpopvtime"), analyzeArgs);
  if (cpopvtimefile_ != 0) {
    if (analyzeArgs.hasKey("normpop"))
      norm_pop_ = Node::CLUSTERPOP;
    else if (analyzeArgs.hasKey("normframe"))
      norm_pop_ = Node::FRAME;
    else
      norm_pop_ = Node::NONE;
  }

  // Cluster number vs time 
  DataFile* cnumvtimefile = DFL.AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  // NOTE overall set name extracted here. This should be the last thing done.
  dsname_ = analyzeArgs.GetStringNext();
  cnumvtime_ = DSL.AddSet(DataSet::INTEGER, dsname_, "Cnum");
  if (cnumvtime_ == 0) return 1;
  if (cnumvtimefile != 0) cnumvtimefile->AddDataSet( cnumvtime_ );

  return 0;
}

// -----------------------------------------------------------------------------
void Cpptraj::Cluster::Control::Info() const {
  if (metric_    != 0) metric_->Info();
  if (algorithm_ != 0) algorithm_->Info();
}

/** Figure out which frames to cluster, cache distances if necessary, then
  * perform clustering. Afterwards sieved frames will be added in if necessary
  * and representative frames will be determined.
  */
int Cpptraj::Cluster::Control::Run() {
  if (metric_ == 0 || algorithm_ == 0) { // TODO check pmatrix_?
    mprinterr("Internal Error: Cluster::Control is not set up.\n");
    return 1;
  }

  // Timers
  timer_run_.Start();
  timer_setup_.Start();
  // Set up the Metric
  if (metric_->Setup()) {
    mprinterr("Error: Metric setup failed.\n");
    return 1;
  }

  // Figure out which frames to cluster
  frameSieve_.Clear();
  int frameSelectErr = 1;
  switch ( frameSelect_ ) {
    case UNSPECIFIED:
      frameSelectErr = frameSieve_.SetFramesToCluster(sieve_, metric_->Ntotal(), sieveSeed_);
      break;
    case FROM_CACHE :
      mprintf("\tClustering frames present in pairwise cache '%s'\n", pmatrix_.Cache().legend());
      frameSelectErr = frameSieve_.SetupFromCache( pmatrix_.Cache(), metric_->Ntotal() );
      break;
    default :
      mprinterr("Internal Error: Cluster::Control::Run(): Unhandled frame selection type.\n");
  }
  if (frameSelectErr != 0) {
    mprinterr("Error: Cluster frame selection failed.\n");
    return 1;
  } 
  if (verbose_ >= 0) {
    if (frameSieve_.FramesToCluster().size() < metric_->Ntotal())
      mprintf("\tClustering %zu of %u points.\n", frameSieve_.FramesToCluster().size(),
              metric_->Ntotal());
    else
      mprintf("\tClustering %u points.\n", metric_->Ntotal());
  }
  Cframes const& framesToCluster = frameSieve_.FramesToCluster();

  timer_setup_.Stop();
  timer_pairwise_.Start();

  // Cache distances if necessary
  if (pmatrix_.CacheDistances( framesToCluster, sieve_ )) return 1;
  if (pmatrix_.HasCache() && verbose_ > 1)
    pmatrix_.Cache().PrintCached();

  timer_pairwise_.Stop();
  timer_cluster_.Start();

  // Cluster
  if (algorithm_->DoClustering(clusters_, framesToCluster, pmatrix_) != 0) {
    mprinterr("Error: Clustering failed.\n");
    return 1;
  }

  timer_cluster_.Stop();

  // ---------------------------------------------
  timer_post_.Start();
  //Timer cluster_post_coords;
  timer_post_renumber_.Start();
  // Update cluster centroids here in case they need to be used to 
  // restore sieved frames
  clusters_.UpdateCentroids( metric_ );

  // Add sieved frames to existing clusters.
  if ( sieveRestore_ != NO_RESTORE ) {
    // Restore sieved frames
    mprintf("\tRestoring sieved frames.\n");
    switch (sieveRestore_) {
      case CLOSEST_CENTROID :
        clusters_.AddFramesByCentroid( frameSieve_.SievedOut(), metric_ ); break;
      case EPSILON_CENTROID :
      case EPSILON_FRAME    :
        clusters_.AddFramesByCentroid( frameSieve_.SievedOut(), metric_,
                                       (sieveRestore_ == EPSILON_CENTROID),
                                       restoreEpsilon_ );
        break;
      default:
        mprinterr("Internal Error: Unhandled sieve restore type.\n");
        return 1;
    }
    // Re-calculate the cluster centroids
    clusters_.UpdateCentroids( metric_ );
  }

  // Sort by population and renumber
  clusters_.Sort();

  timer_post_renumber_.Stop();
  timer_post_bestrep_.Start();

  // Find best representative frames for each cluster.
  if (BestReps::FindBestRepFrames(bestRep_, nRepsToSave_, clusters_, pmatrix_,
                                  frameSieve_.SievedOut(), verbose_))
  {
    mprinterr("Error: Finding best representative frames for clusters failed.\n");
    return 1;
  }

  timer_post_bestrep_.Stop();

  // DEBUG - print clusters to stdout
  if (verbose_ > 0) {
    mprintf("\nFINAL CLUSTERS:\n");
    clusters_.PrintClusters();
  }

  // TODO assign reference names
  timer_run_.Stop();
  return 0;
}

int Cpptraj::Cluster::Control::Output(DataSetList& DSL) {
  timer_output_.Start();
  // Info
  if (!suppressInfo_) {
    CpptrajFile outfile;
    if (outfile.OpenWrite( clusterinfo_ )) return 1;
    timer_output_info_.Start();
    Output::PrintClustersToFile(outfile, clusters_, *algorithm_, metric_, 
                                frameSieve_.SieveValue(), frameSieve_.FramesToCluster());
    timer_output_info_.Stop();
    outfile.CloseFile();
  }

  // Silhouette
  if (!sil_file_.empty()) {
    if (frameSieve_.SieveValue() != 1 && !includeSieveInCalc_)
      mprintf("Warning: Silhouettes do not include sieved frames.\n");
    clusters_.CalcSilhouette(pmatrix_, frameSieve_.SievedOut(), includeSieveInCalc_);
    CpptrajFile Ffile, Cfile;
    if (Ffile.OpenWrite(sil_file_ + ".frame.dat")) return 1;
    Output::PrintSilhouetteFrames(Ffile, clusters_);
    Ffile.CloseFile();
    if (Cfile.OpenWrite(sil_file_ + ".cluster.dat")) return 1;
    Output::PrintSilhouettes(Cfile, clusters_);
    Cfile.CloseFile();
  }

  // Print a summary of clusters
  if (!summaryfile_.empty()) {
    timer_output_summary_.Start();
    CpptrajFile outfile;
    if (outfile.OpenWrite(summaryfile_)) {
      mprinterr("Error: Could not set up cluster summary file.\n");
      return 1;
    }
    Output::Summary(outfile, clusters_, *algorithm_, pmatrix_, includeSieveInCalc_,
                    frameSieve_.SievedOut());
    timer_output_summary_.Stop();
  }

  // Cluster number vs time
  if (cnumvtime_ != 0) {
    if (clusters_.CreateCnumVsTime( (DataSet_integer*)cnumvtime_, metric_->Ntotal() ))
    {
      mprinterr("Error: Creation of cluster num vs time data set failed.\n");
      return 1;
    }
  }

  // Cluster population vs time
  if (cpopvtimefile_ != 0) {
    MetaData md( dsname_, "Pop" );
    DataSet::SizeArray setsize(1, metric_->Ntotal());
    for (List::cluster_iterator node = clusters_.begin();
                                node != clusters_.end(); ++node)
    {
      md.SetIdx( node->Num() );
      DataSet_float* ds = (DataSet_float*)DSL.AddSet( DataSet::FLOAT, md );
      if (ds == 0) {
        mprinterr("Error: Could not allocate cluster pop vs time DataSet\n");
        return 1;
      }
      ds->Allocate( setsize );
      if (cpopvtimefile_ != 0) cpopvtimefile_->AddDataSet( ds );
      node->CalcCpopVsTime( *ds, metric_->Ntotal(), norm_pop_ );
    }
  }

  timer_output_.Stop();
  return 0;
}

void Cpptraj::Cluster::Control::Timing(double ttotal) const {
  mprintf("\tCluster timing data:\n");
  // Run Timing data
  timer_setup_.WriteTiming(       2, "Cluster Init.  :", timer_run_.Total());
  timer_pairwise_.WriteTiming(    2, "Pairwise Calc. :", timer_run_.Total());
  timer_cluster_.WriteTiming(     2, "Clustering     :", timer_run_.Total());
  algorithm_->Timing( timer_cluster_.Total() );
  timer_post_.WriteTiming(        2, "Cluster Post.  :", timer_run_.Total());
  timer_post_renumber_.WriteTiming( 3, "Cluster renumbering/sieve restore :", timer_post_.Total());
  timer_post_bestrep_.WriteTiming(  3, "Find best rep.                    :", timer_post_.Total());
  timer_run_.WriteTiming(   1, "Run Total    :", ttotal);
  // Output Timing data
  timer_output_info_.WriteTiming(   2, "Info calc      :", timer_output_.Total());
  timer_output_summary_.WriteTiming(2, "Summary calc   :", timer_output_.Total());
  //cluster_post_coords.WriteTiming(2, "Coordinate writes", cluster_post.Total());
  timer_output_.WriteTiming(1, "Output Total :", ttotal);
}

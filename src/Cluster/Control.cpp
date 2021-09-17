#include "Control.h"
#include "Output.h"
#include "../ArgList.h"
#include "../BufferedLine.h" // For loading info file
#include "../CpptrajStdio.h"
#include "../DataFile.h" // For loading pairwise cache
#include "../DataFileList.h"
#include "../DataSet_Coords.h"
#include "../DataSet_float.h"
#include "../DataSet_integer.h"
#include "../DataSet_PairwiseCache.h"
// Algorithms
#include "Algorithm_HierAgglo.h"
#include "Algorithm_DBscan.h"
#include "Algorithm_Kmeans.h"
#include "Algorithm_DPeaks.h"
// Results
#include "Results_Coords.h"

/** CONSTRUCTOR */
Cpptraj::Cluster::Control::Control() :
  algorithm_(0),
  results_(0),
  verbose_(0),
  frameSelect_(UNSPECIFIED),
  sieve_(1),
  sieveSeed_(-1),
  sieveRestore_(NO_RESTORE),
  restoreEpsilon_(0.0),
  includeSieveInCalc_(false),
  includeSieveCdist_(false),
  bestRep_(BestReps::NO_REPS),
  nRepsToSave_(1),
  suppressInfo_(false),
  cnumvtime_(0),
  grace_color_(false),
  clustersVtime_(0),
  windowSize_(0),
  cpopvtimefile_(0),
  norm_pop_(Node::NONE),
  calc_lifetimes_(false),
  drawGraph_(NO_DRAWGRAPH),
  draw_tol_(0),
  draw_maxit_(0),
  debug_(0),
  metricContribFile_(0)
{}

/** DESTRUCTOR */
Cpptraj::Cluster::Control::~Control() {
  if (algorithm_ != 0) delete algorithm_;
  if (results_ != 0  ) delete results_;
}

// -----------------------------------------------------------------------------
/** \return Pointer to Algorithm of given type. */
Cpptraj::Cluster::Algorithm* Cpptraj::Cluster::Control::AllocateAlgorithm(Algorithm::AType atype)
{
  Algorithm* alg = 0;
  switch (atype) {
    case Algorithm::HIERAGGLO : alg = new Algorithm_HierAgglo(); break;
    case Algorithm::DBSCAN    : alg = new Algorithm_DBscan(); break;
    case Algorithm::KMEANS    : alg = new Algorithm_Kmeans(); break;
    case Algorithm::DPEAKS    : alg = new Algorithm_DPeaks(); break;
    default : mprinterr("Error: Unhandled Algorithm in AllocateAlgorithm.\n");
  }
  return alg;
}

/** Recognized algorithm keywords. */
const char* Cpptraj::Cluster::Control::AlgorithmArgs_ =
  "{hieragglo|dbscan|kmeans|dpeaks}";

/** Set up Algorithm from keyword + arguments. */
int Cpptraj::Cluster::Control::AllocateAlgorithm(ArgList& analyzeArgs) {
  Algorithm::AType atype = Algorithm::UNSPECIFIED;
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
static int Err(int code) {
  switch (code) {
    case 0: mprinterr("Error: Could not open info file.\n"); break;
    case 1: mprinterr("Error: Unexpected end of info file.\n"); break;
    case 2: mprinterr("Error: Invalid number of clusters in info file.\n"); break;
    case 3: mprinterr("Error: Invalid number of frames in info file.\n"); break;
  }
  return 1;
}

/** Read clustering info from existing info file. */
int Cpptraj::Cluster::Control::ReadInfo(std::string const& fname) {
  if (fname.empty()) {
    mprinterr("Error: No cluster info filename given.\n");
    return 1;
  }
  mprintf("\tReading clusters from info file '%s'\n", fname.c_str());
  BufferedLine infile;
  if (infile.OpenFileRead( fname )) return Err(0);
  const char* ptr = infile.Line();
  if (ptr == 0) return Err(1);
  ArgList infoLine( ptr, " " );
  int nclusters = infoLine.getKeyInt("#Clustering:", -1);
  if (nclusters == -1) return Err(2);
  int nframes = infoLine.getKeyInt("clusters", -1);
  if (nframes == -1) return Err(3);
  mprintf("\tNumber of frames in info file: %i\n", nframes);
//  if (nframes != (int)FrameDistances().OriginalNframes()) {
//    mprinterr("Error: # frames in cluster info file (%i) does not match"
//              " current # frames (%zu)\n", nframes, FrameDistances().OriginalNframes());
//    return 1;
//  }
  // Scan down to clusters
  std::string algorithmStr;
  while (ptr[0] == '#') {
    ptr = infile.Line();
    if (ptr == 0) return Err(1);
    // Save previous clustering info. Includes newline.
    if (ptr[1] == 'A' && ptr[2] == 'l' && ptr[3] == 'g')
      algorithmStr.assign( ptr + 12 ); // Right past '#Algorithm: '
  }
  if (!algorithmStr.empty()) mprintf("\tAlgorithm in info file: %s\n", algorithmStr.c_str());
  // Read clusters
  Cframes frames;
  for (int cnum = 0; cnum != nclusters; cnum++) {
    if (ptr == 0) return Err(1);
    frames.clear();
    // TODO: Check for busted lines?
    for (int fidx = 0; fidx != nframes; fidx++) {
      if (ptr[fidx] == 'X')
        frames.push_back( fidx );
    }
    clusters_.AddCluster( Node(metrics_, frames, cnum) );
    mprintf("\tRead cluster %i, %zu frames.\n", cnum, frames.size());
    ptr = infile.Line();
  }
  infile.CloseFile();
  return 0;
}

/** Set up clusters from a cluster number vs time data set. */
int Cpptraj::Cluster::Control::InitClustersFromSet(DataSet* cnvt) {
  if (cnvt == 0) {
    mprinterr("Internal Error: InitClustersFromSet: Null set.\n");
    return 1;
  }
  mprintf("\tReading clusters from cluster num vs time set '%s'\n", cnvt->legend());
  if (cnvt->Group() != DataSet::SCALAR_1D) {
    mprinterr("Error: Set '%s' is not a scalar 1D set.\n", cnvt->legend());
    return 1;
  }
  if (cnvt->Type() != DataSet::INTEGER) {
    mprintf("Warning: Set '%s' is not an integer set. Floating point values\n"
            "Warning:  will be rounded to integer values.\n", cnvt->legend());
  }
  // Hold frames for each cluster
  std::vector<Cframes> clusterFrames;
  DataSet_1D const& ds = static_cast<DataSet_1D const&>( *cnvt );
  for (unsigned int idx = 0; idx != cnvt->Size(); idx++)
  {
    int cnum = ds.Dval(idx);
    if (cnum > -1) {
      if (cnum >= (int)clusterFrames.size()) {
        while ((int)clusterFrames.size() <= cnum)
          clusterFrames.push_back( Cframes() );
      }
      clusterFrames[cnum].push_back( (int)idx );
    }
  }
  int clusterNum = 0;
  for (std::vector<Cframes>::const_iterator frames = clusterFrames.begin();
                                            frames != clusterFrames.end();
                                          ++frames, ++clusterNum)
    clusters_.AddCluster( Node(metrics_, *frames, clusterNum) );
  return 0;
}

// -----------------------------------------------------------------------------

const char* Cpptraj::Cluster::Control::SieveArgs1_ =
  "[sieve <#> [sieveseed <#>] [random] [includesieveincalc] [includesieved_cdist]";

const char* Cpptraj::Cluster::Control::SieveArgs2_ =
  " [{sievetoframe|sievetocentroid|closestcentroid}] [repsilon <restore epsilon>]]";

const char* Cpptraj::Cluster::Control::BestRepArgs_ =
  "[bestrep {cumulative|centroid|cumulative_nosieve}] [savenreps <#>]";

const char* Cpptraj::Cluster::Control::OutputArgs1_ =
  "[out <cnumvtime> [gracecolor]] [noinfo|info <file>] [summary <file>]";

const char* Cpptraj::Cluster::Control::OutputArgs2_ =
  "[summarysplit <splitfile>] [splitframe <comma-separated frame list>]";

const char* Cpptraj::Cluster::Control::OutputArgs3_ =
  "[clustersvtime <file> [cvtwindow <#>]] [sil <prefix>] [metricstats <file>]";

const char* Cpptraj::Cluster::Control::OutputArgs4_ =
  "[cpopvtime <file> [{normpop|normframe}]] [lifetime]";

const char* Cpptraj::Cluster::Control::GraphArgs_ =
  "[{drawgraph|drawgraph3d} [draw_tol <tolerance>] [draw_maxit <iterations]]";

/** Common setup. */
//int Cpptraj::Cluster::Control::Common(ArgList& analyzeArgs, DataSetList& DSL, DataFileList& DFL)
int Cpptraj::Cluster::Control::SetupClustering(DataSetList const& setsToCluster,
                                               DataSet* coordsSet,
                                               ArgList& analyzeArgs,
                                               DataSetList& DSL,
                                               DataFileList& DFL,
                                               int verboseIn)
{
  verbose_ = verboseIn;
  clusters_.SetDebug( verbose_ );

  if (setsToCluster.empty()) {
    mprinterr("Error: No sets to cluster.\n");
    return 1;
  }

  // Initialize metrics based on set types. Also initializes the pairwise cache.
  if (metrics_.Initialize( setsToCluster, DSL, DFL, analyzeArgs, verbose_ )) {
    mprinterr("Error: Could not initialize distance metrics/cache.\n");
    return 1;
  }

  // Determine if metric contributions should be determined
  metricContribFile_ = 0;
  std::string metricContribArg = analyzeArgs.GetStringKey("metricstats");
  if (!metricContribArg.empty()) {
    if (setsToCluster.size() < 2) {
      mprintf("Warning: 'metricstats' is only relevant with 2 or more input data sets. Skipping.\n");
    } else {
      metricContribFile_ = DFL.AddCpptrajFile(metricContribArg, "Metric Contributions",
                                              DataFileList::TEXT);
      if (metricContribFile_ == 0) {
        mprinterr("Error: Could not set up 'metricstats' file.\n");
        return 1;
      }
    }
  }

  // Initialize clusters from existing info file. Metric must already be set up.
  if (analyzeArgs.hasKey("readinfo") ||
      analyzeArgs.hasKey("readtxt"))
  {
    std::string iname = analyzeArgs.GetStringKey("infofile");
    if (!iname.empty()) {
      if (ReadInfo( iname )) return 1;
    } else {
      iname = analyzeArgs.GetStringKey("cnvtset");
      if (!iname.empty()) {
        DataSet* cnvt = DSL.GetDataSet( iname );
        if (cnvt == 0) {
          mprinterr("Error: Input cluster num vs time set '%s' not found.\n", iname.c_str());
          return 1;
        }
        if (InitClustersFromSet( cnvt )) return 1;
      } else {
        mprinterr("Error: Must specify either 'infofile' or 'cnvtset' with 'readinfo'\n");
        return 1;
      }
    }
  }

  // Set up results that depend on COORDS DataSet
  if (coordsSet != 0) {
    if (results_ != 0) delete results_;
    results_ = (Results*)new Results_Coords( static_cast<DataSet_Coords*>( coordsSet ) );
    if (results_ == 0) return 1;
    if (results_->GetOptions(analyzeArgs, DSL, metrics_)) return 1;
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
  includeSieveCdist_ = analyzeArgs.hasKey("includesieved_cdist");
  if (includeSieveCdist_)
    mprintf("Warning: 'includesieved_cdist' may be very slow.\n");

  // Determine how frames to cluster will be chosen
  if (frameSelect_ == UNSPECIFIED) {
    // If no other frame selection option like sieve provided and an already
    // set up cache is present, use the cached frames.
    if (sieve_ == 1 && metrics_.HasCache() && metrics_.Cache().Size() > 0)
    {
      frameSelect_ = FROM_CACHE;
      if (metrics_.Cache().SieveVal() != 1) {
        sieve_ = metrics_.Cache().SieveVal();
        mprintf("Warning: No sieve specified; using sieve value from cache '%s': %i\n",
                metrics_.Cache().legend(), sieve_);
      }
    }
  }

  // Choose sieve restore options
  if (sieve_ != 1)
  {
    // Determine sieve restore type
    if (sieveRestore_ == NO_RESTORE) {
      // No option set yet. See if a keyword has been specified.
      if (analyzeArgs.hasKey("nosieverestore")) // Hidden option, currently for testing only
        sieveRestore_ = NO_RESTORE;
      else if (analyzeArgs.hasKey("sievetoframe"))
        sieveRestore_ = EPSILON_FRAME;
      else if (analyzeArgs.hasKey("sievetocentroid"))
        sieveRestore_ = EPSILON_CENTROID;
      else if (analyzeArgs.hasKey("closestcentroid"))
        sieveRestore_ = CLOSEST_CENTROID;
      else {
        // Nothing specified yet. Choose restore option based on algorithm.
        if (algorithm_->Type() == Algorithm::DBSCAN ||
            algorithm_->Type() == Algorithm::DPEAKS)
          sieveRestore_ = EPSILON_CENTROID;
        else
          sieveRestore_ = CLOSEST_CENTROID;
      }
    }
    // Determine sieve restore epsilon
    if (sieveRestore_ == EPSILON_FRAME ||
        sieveRestore_ == EPSILON_CENTROID)
    {
      // Epsilon-based sieve restore
      if (algorithm_->Type() == Algorithm::DBSCAN ||
          algorithm_->Type() == Algorithm::DPEAKS)
      {
        // Using a density-based algorithm with epsilon-based restore;
        // use restore epsilon from algorithm.
        restoreEpsilon_ = algorithm_->Epsilon();
      }
    }
    double rEps = analyzeArgs.getKeyDouble("repsilon", -1.0);
    if (rEps > 0)
      restoreEpsilon_ = rEps;
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

  // Draw graph
  if (analyzeArgs.hasKey("drawgraph"))
    drawGraph_ = TWOD;
  else if (analyzeArgs.hasKey("drawgraph3d"))
    drawGraph_ = THREED;
  else
    drawGraph_ = NO_DRAWGRAPH;
  draw_maxit_ = analyzeArgs.getKeyInt("draw_maxit", 1000);
  draw_tol_ = analyzeArgs.getKeyDouble("draw_tol", 1.0E-5);

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

  // Number of unique clusters vs time
  DataFile* clustersvtimefile = DFL.AddDataFile(analyzeArgs.GetStringKey("clustersvtime"),
                                                analyzeArgs);
  windowSize_ = analyzeArgs.getKeyInt("cvtwindow", 0);

  // Create cluster lifetime data sets?
  calc_lifetimes_ = analyzeArgs.hasKey("lifetime");

  // Cluster number vs time
  grace_color_ = analyzeArgs.hasKey("gracecolor"); 
  DataFile* cnumvtimefile = DFL.AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);

  // Cluster split analysis
  splitfile_ = analyzeArgs.GetStringKey("summarysplit");
  if (splitfile_.empty()) // For backwards compat.
    splitfile_ = analyzeArgs.GetStringKey("summaryhalf");
  if (!splitfile_.empty()) {
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
        return 1;
      }
    }
  }

  // Overall set name extracted here. All other arguments should already be processed. 
  dsname_ = analyzeArgs.GetStringNext();
  if (dsname_.empty())
    dsname_ = DSL.GenerateDefaultName("CLUSTER");
  // ---------------------------------------------

  // Cluster number vs time data set
  cnumvtime_ = DSL.AddSet(DataSet::INTEGER, dsname_, "Cnum");
  if (cnumvtime_ == 0) return 1;
  if (cnumvtimefile != 0) cnumvtimefile->AddDataSet( cnumvtime_ );

  // If cache was allocated, add to the DataSetList so it is after cnumvtime for pytraj
  if (metrics_.CacheWasAllocated())
    DSL.AddSet( metrics_.CachePtr() );

  // Set up number of unique clusters vs time DataSet
  if (clustersvtimefile != 0) {
    if (windowSize_ < 2) {
      mprinterr("Error: For # clusters seen vs time, cvtwindow must be specified and > 1\n");
      return 1;
    }
    clustersVtime_ = DSL.AddSet(DataSet::INTEGER, MetaData(cnumvtime_->Meta().Name(), "NCVT"));
    if (clustersVtime_ == 0) return 1;
    clustersvtimefile->AddDataSet( clustersVtime_ );
  }

  return 0;
}

/** Print help text to STDOUT. */
void Cpptraj::Cluster::Control::Help() {
  mprintf("\t[<name>] [<Algorithm>] [<Metric>] [<Pairwise>] [<Sieve>] [<BestRep>]\n"
          "\t[<Output>] [<Coord. Output>] [<Graph>]\n"
          "\t[readinfo {infofile <info file> | cnvtset <dataset>}]\n");
  mprintf("  Algorithm Args: [%s]\n", AlgorithmArgs_);
  Algorithm_HierAgglo::Help();
  Algorithm_DBscan::Help();
  Algorithm_Kmeans::Help();
  Algorithm_DPeaks::Help();
  mprintf("  Metric Args:\n");
  mprintf("\t%s\n", MetricArray::MetricArgs_);
  mprintf("  Pairwise Args:\n");
  mprintf("\t%s\n", MetricArray::PairwiseArgs1_);
  mprintf("\t%s\n", MetricArray::PairwiseArgs2_);
  mprintf("  Sieve Args:\n");
  mprintf("\t%s\n", SieveArgs1_);
  mprintf("\t%s\n", SieveArgs2_);
  mprintf("  BestRep Args:\n");
  mprintf("\t%s\n", BestRepArgs_);
  mprintf("  Output Args:\n");
  mprintf("\t%s\n", OutputArgs1_);
  mprintf("\t%s\n", OutputArgs2_);
  mprintf("\t%s\n", OutputArgs3_);
  mprintf("\t%s\n", OutputArgs4_);
  mprintf("  Coordinate Output Args:\n");
  Results_Coords::Help();
  mprintf("  Graph Args:\n");
  mprintf("\t%s\n", GraphArgs_);
  mprintf("  Cluster the data from given input data set(s) (currently any combination\n"
          "  of COORDS and 1D sets are supported) using the specified clustering\n"
          "  algorithm. If any input sets are COORDS sets, the metric used to calculate\n"
          "  frame-frame distances can also be specified (default RMSD).\n"
          "  In order to speed up clustering of large sets, the 'sieve' keyword can be\n"
          "  used to cluster a subset of the data, after which the remaining frames\n"
          "  are added back in. In addition, subsequent clustering calculations can be\n"
          "  sped up by writing to/reading from a pair-wise distance cache file/DataSet\n"
          "  specified by 'pairdist'.\n");
}

// -----------------------------------------------------------------------------
void Cpptraj::Cluster::Control::Info() const {
  metrics_.Info();
  if (metricContribFile_ != 0)
    mprintf("\tContributions of each metric to the total distance will be printed to '%s'.\n",
            metricContribFile_->Filename().full());

  if (algorithm_ != 0) algorithm_->Info();

  if (results_   != 0)
    results_->Info();
  else
    mprintf("\tNo coordinates set provided for cluster results.\n");

  if (frameSelect_ == FROM_CACHE) {
    mprintf("\tWill cluster frames present in existing pairwise cache '%s'\n",
            metrics_.Cache().legend());
  } else {
    if (sieve_ > 1)
      mprintf("\tInitial clustering sieve value is %i frames.\n", sieve_);
    else if (sieve_ < -1) {
      mprintf("\tInitial clustering will be randomly sieved (with value %i)", -sieve_);
      if (sieveSeed_ > 0) mprintf(" using random seed %i", sieveSeed_);
      mprintf(".\n");
    }
  }

  if (sieveRestore_ != NO_RESTORE) {
    mprintf("\tRestoring sieved frames");
    if (sieveRestore_ == CLOSEST_CENTROID)
      mprintf(" by closest distance to centroid.\n");
    else if (sieveRestore_ == EPSILON_CENTROID)
      mprintf(" if within epsilon %f of a centroid (less accurate but faster).\n", restoreEpsilon_);
    else if (sieveRestore_ == EPSILON_FRAME)
      mprintf(" if within epsilon %f of a frame (more accurate and identifies noise but slower).\n", restoreEpsilon_);
  }

  mprintf("\tRepresentative frames will be chosen by");
  switch (bestRep_) {
    case BestReps::CUMULATIVE: mprintf(" lowest cumulative distance to all other frames.\n"); break;
    case BestReps::CENTROID  : mprintf(" closest distance to cluster centroid.\n"); break;
    case BestReps::CUMULATIVE_NOSIEVE:
      mprintf(" lowest cumulative distance to all other frames (ignore sieved frames).\n");
      break;
    default: // sanity check
      mprintf("\n");
  }
  if (nRepsToSave_ > 1)
    mprintf("\tThe top %i representative frames will be determined.\n", nRepsToSave_);

  if (!clusterinfo_.empty())
    mprintf("\tCluster information will be written to %s\n",clusterinfo_.c_str());
  if (!summaryfile_.empty()) {
    mprintf("\tSummary of cluster results will be written to %s\n",summaryfile_.c_str());
    if (sieve_ != 1) {
      if (includeSieveInCalc_)
        mprintf("\tInternal cluster averages will include sieved frames.\n");
      if (includeSieveCdist_)
        mprintf("\tBetween-cluster distances will include sieved frames.\n");
    }
  }
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

  if (cnumvtime_ != 0) {
    mprintf("\tCluster number vs time data set: %s\n", cnumvtime_->legend());
    if (grace_color_)
      mprintf("\tGrace color instead of cluster number (1-15, 0 = noise) will be saved.\n");
  }

  if (clustersVtime_ != 0) {
    mprintf("\tNumber of unique clusters observed over windows of size %i data set: %s\n",
            windowSize_, clustersVtime_->legend());
  }

  if (cpopvtimefile_ != 0) {
    mprintf("\tCluster pop vs time will be written to %s", cpopvtimefile_->DataFilename().base());
    if (norm_pop_ == Node::CLUSTERPOP)
      mprintf(" (normalized by cluster size)");
    else if (norm_pop_ == Node::FRAME)
      mprintf(" (normalized by frame)");
    mprintf("\n");
  }

  if (calc_lifetimes_)
    mprintf("\tCluster lifetime data sets will be created with aspect 'Lifetime'\n");

  if (!splitfile_.empty()) {
    mprintf("\tSummary comparing parts of trajectory data for clusters will be written to %s\n",
            splitfile_.c_str());
    if (!splitFrames_.empty()) {
      mprintf("\t\tFrames will be split at:");
      for (std::vector<int>::const_iterator f = splitFrames_.begin(); f != splitFrames_.end(); ++f)
        mprintf(" %i", *f);
      mprintf("\n");
    } else
      mprintf("\t\tFrames will be split at the halfway point.\n");
  }
  // TODO loaded refs_
  if (drawGraph_ != NO_DRAWGRAPH)
    mprintf("\tEXPERIMENTAL: Force-directed graph will be drawn from pairwise distances.\n"
            "\t              Max iterations= %i, min tolerance= %g\n",
                             draw_maxit_, draw_tol_);
}

/** Figure out which frames to cluster, cache distances if necessary, then
  * perform clustering. Afterwards sieved frames will be added in if necessary
  * and representative frames will be determined.
  */
int Cpptraj::Cluster::Control::Run() {
  if (metrics_.empty() || algorithm_ == 0) {
    mprinterr("Internal Error: Cluster::Control is not set up.\n");
    return 1;
  }

  // Timers
  timer_run_.Start();
  timer_setup_.Start();
  // Set up the Metric
  if (metrics_.Setup()) {
    mprinterr("Error: Metric setup failed.\n");
    return 1;
  }
  if (metrics_.Ntotal() < 1) {
    mprinterr("Error: No frames to cluster; all input sets are empty.\n");
    return 1;
  }

  // Figure out which frames to cluster
  frameSieve_.Clear();
  int frameSelectErr = 1;
  switch ( frameSelect_ ) {
    case UNSPECIFIED:
      frameSelectErr = frameSieve_.SetFramesToCluster(sieve_, metrics_.Ntotal(), sieveSeed_);
      break;
    case FROM_CACHE :
      mprintf("\tClustering frames present in pairwise cache '%s'\n", metrics_.Cache().legend());
      frameSelectErr = frameSieve_.SetupFromCache( metrics_.Cache(), metrics_.Ntotal() );
      break;
    default :
      mprinterr("Internal Error: Cluster::Control::Run(): Unhandled frame selection type.\n");
  }
  if (frameSelectErr != 0) {
    mprinterr("Error: Cluster frame selection failed.\n");
    return 1;
  } 
  if (verbose_ >= 0) {
    if (frameSieve_.FramesToCluster().size() < metrics_.Ntotal())
      mprintf("\tClustering %zu of %u points.\n", frameSieve_.FramesToCluster().size(),
              metrics_.Ntotal());
    else
      mprintf("\tClustering %u points.\n", metrics_.Ntotal());
  }
  Cframes const& framesToCluster = frameSieve_.FramesToCluster();

  // Calculate contribution to total for each metric
  if (metricContribFile_ != 0)
    metrics_.CalculateMetricContributions(framesToCluster, *metricContribFile_);

  timer_setup_.Stop();
  timer_pairwise_.Start();

  // Cache distances if necessary
  if (metrics_.CacheDistances( framesToCluster, sieve_ )) return 1;
  if (metrics_.HasCache() && verbose_ > 1)
    metrics_.Cache().PrintCached();

  timer_pairwise_.Stop();
  timer_cluster_.Start();

  // Cluster
  if (algorithm_->DoClustering(clusters_, framesToCluster, metrics_) != 0) {
    mprinterr("Error: Clustering failed.\n");
    return 1;
  }

  timer_cluster_.Stop();

  // ---------------------------------------------
  if (clusters_.Nclusters() == 0) {
    mprintf("\tNo clusters found.\n");
  } else {
    timer_post_.Start();
    //Timer cluster_post_coords;
    timer_post_renumber_.Start();
    // Update cluster centroids here in case they need to be used to 
    // restore sieved frames
    clusters_.UpdateCentroids( metrics_ );

    // Add sieved frames to existing clusters.
    if ( sieveRestore_ != NO_RESTORE ) {
      // Restore sieved frames
      mprintf("\tRestoring sieved frames.\n");
      switch (sieveRestore_) {
        case CLOSEST_CENTROID :
          clusters_.AddFramesByCentroid( frameSieve_.SievedOut(), metrics_ ); break;
        case EPSILON_CENTROID :
        case EPSILON_FRAME    :
          clusters_.AddFramesByCentroid( frameSieve_.SievedOut(), metrics_,
                                         (sieveRestore_ == EPSILON_CENTROID),
                                         restoreEpsilon_ );
          break;
        default:
          mprinterr("Internal Error: Unhandled sieve restore type.\n");
          return 1;
      }
      // Re-calculate the cluster centroids
      clusters_.UpdateCentroids( metrics_ );
    }

    // Sort by population and renumber
    clusters_.Sort();

    timer_post_renumber_.Stop();
    timer_post_bestrep_.Start();

    // Find best representative frames for each cluster.
    BestReps findBestReps;
    if (findBestReps.InitBestReps(bestRep_, nRepsToSave_, verbose_)) {
      mprinterr("Error: Initializing best representative frames search failed.\n");
      return 1;
    }
    if (findBestReps.FindBestRepFrames(clusters_, metrics_, frameSieve_.SievedOut())) {
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
  }
  timer_run_.Stop();
  return 0;
}

/** Write results to files etc. */
int Cpptraj::Cluster::Control::Output(DataSetList& DSL) {
  if (clusters_.Nclusters() == 0) return 0;
  timer_output_.Start();

  // Results that require additional calculations
  if (results_ != 0) {
    timer_output_results_.Start();
    results_->CalcResults( clusters_ );
    timer_output_results_.Stop();
  }

  // Info
  if (!suppressInfo_) {
    CpptrajFile outfile;
    if (outfile.OpenWrite( clusterinfo_ )) return 1;
    timer_output_info_.Start();
    Output::PrintClustersToFile(outfile, clusters_, *algorithm_, metrics_, 
                                frameSieve_.SieveValue(), frameSieve_.FramesToCluster());
    timer_output_info_.Stop();
    outfile.CloseFile();
  }

  // Silhouette
  if (!sil_file_.empty()) {
    if (frameSieve_.SieveValue() != 1 && !includeSieveInCalc_)
      mprintf("Warning: Silhouettes do not include sieved frames.\n");
    clusters_.CalcSilhouette(metrics_, frameSieve_.SievedOut(), includeSieveInCalc_);
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
    Output::Summary(outfile, clusters_, *algorithm_, metrics_, includeSieveInCalc_,
                    includeSieveCdist_, frameSieve_.SievedOut());
    timer_output_summary_.Stop();
  }

  // Summary by parts
  if (!splitfile_.empty()) {
    CpptrajFile outfile;
    if (outfile.OpenWrite(splitfile_)) {
      mprinterr("Error: Could not set up summary split file.\n");
      return 1;
    }
    BestReps findBestReps;
    // TODO support > 1 best rep?
    if (findBestReps.InitBestReps(bestRep_, 1, verbose_)) {
      mprinterr("Error: Init of best reps calc for summary by parts failed.\n");
      return 1;
    }
    // TODO just pass in metrics_
    Output::Summary_Part(outfile, metrics_.Ntotal(), splitFrames_, clusters_,
                         findBestReps, metrics_, frameSieve_.SievedOut());
  }

  // Cluster number vs time
  if (cnumvtime_ != 0) {
    int err = 0;
    if (grace_color_)
      err = clusters_.CreateCnumVsTime(*((DataSet_integer*)cnumvtime_), metrics_.Ntotal(), 1, 15);
    else
      err = clusters_.CreateCnumVsTime(*((DataSet_integer*)cnumvtime_), metrics_.Ntotal(), 0, -1);
    if (err != 0) {
      mprinterr("Error: Creation of cluster num vs time data set failed.\n");
      return 1;
    }
  }

  // Draw cluster Graph
  if (drawGraph_ != NO_DRAWGRAPH) {
    DrawGraph( frameSieve_.FramesToCluster(), metrics_, drawGraph_, cnumvtime_, draw_tol_, draw_maxit_, debug_ );
  }

  // # unique clusters vs time
  if (clustersVtime_ != 0) {
    if (clusters_.NclustersObserved(*((DataSet_integer*)clustersVtime_), metrics_.Ntotal(),
                                    windowSize_))
    {
      mprinterr("Error: Creation of # unique clusters vs time data set failed.\n");
      return 1;
    }
  }

  // Cluster population vs time
  if (cpopvtimefile_ != 0) {
    MetaData md( dsname_, "Pop" );
    DataSet::SizeArray setsize(1, metrics_.Ntotal());
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
      node->CalcCpopVsTime( *ds, metrics_.Ntotal(), norm_pop_ );
    }
  }

  // Cluster lifetime sets
  if (calc_lifetimes_) {
    MetaData md( dsname_, "Lifetime" );
    DataSet::SizeArray setsize(1, metrics_.Ntotal());
    for (List::cluster_iterator node = clusters_.begin();
                                node != clusters_.end(); ++node)
    {
      md.SetIdx( node->Num() );
      DataSet_integer* ds = (DataSet_integer*)DSL.AddSet( DataSet::INTEGER, md );
      if (ds == 0) {
        mprinterr("Error: Could not allocate cluster lifetime DataSet\n");
        return 1;
      }
      ds->Allocate( setsize );
      node->CreateLifetimeSet( *ds, metrics_.Ntotal() );
    }
  }

  // Any other results
  if (results_ != 0) {
    timer_output_results_.Start();
    results_->DoOutput( clusters_ );
    timer_output_results_.Stop();
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
  timer_output_results_.WriteTiming(2, "Results Output :", timer_output_.Total());
  timer_output_.WriteTiming(1, "Output Total :", ttotal);
}

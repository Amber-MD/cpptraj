#include "Control.h"
#include "../CpptrajStdio.h"
#include "../DataSet_Coords.h"
#include "Sieve.h"
// PairwiseMatrix classes
#include "PairwiseMatrix_MEM.h"
// Metric classes
#include "Metric_RMS.h"
// Algorithms
#include "Algorithm_HierAgglo.h"
#include "Algorithm_DBscan.h"

Cpptraj::Cluster::Control::Control() :
  metric_(0),
  pmatrix_(0),
  algorithm_(0),
  verbose_(0),
  sieve_(1),
  sieveSeed_(-1),
  sieveRestore_(NO_RESTORE),
  restoreEpsilon_(0.0),
  bestRep_(BestReps::NO_REPS),
  nRepsToSave_(1)
{}

// -----------------------------------------------------------------------------
/** \return pointer to PairwiseMatrix of specified type. */
Cpptraj::Cluster::PairwiseMatrix* 
  Cpptraj::Cluster::Control::AllocatePairwise(PairwiseMatrix::Type ptype, Metric* metric)
{
  PairwiseMatrix* pmatrix = 0;
  // TODO check for null metric
  switch (ptype) {
    case PairwiseMatrix::MEM : pmatrix = new PairwiseMatrix_MEM(metric); break;
    default: mprinterr("Error: Unhandled PairwiseMatrix in AllocatePairwise.\n");
  }
  return pmatrix;
}

const char* Cpptraj::Cluster::Control::PairwiseArgs =
  "pairwisecache {mem|disk|none}";

/** Set up PairwiseMatrix from arguments. */
int Cpptraj::Cluster::Control::AllocatePairwise(ArgList& analyzeArgs, Metric* metricIn)
{
  if (metricIn == 0) return 1;
  PairwiseMatrix::Type pw_type = PairwiseMatrix::MEM;
  std::string pw_typeString = analyzeArgs.GetStringKey("pairwisecache");
  if (!pw_typeString.empty()) {
    if (pw_typeString == "mem")
      pw_type = PairwiseMatrix::MEM;
    else if (pw_typeString == "disk")
      pw_type = PairwiseMatrix::DISK;
    else if (pw_typeString == "none")
      pw_type = PairwiseMatrix::NOCACHE;
    else {
      mprinterr("Error: Unrecognized option for 'pairwisecache' ('%s')\n", pw_typeString.c_str());
      return 1;
    }
  }
  if (pmatrix_ != 0) delete pmatrix_;
  pmatrix_ = AllocatePairwise( pw_type, metricIn );
  if (pmatrix_ == 0) return 1;
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

const char* Cpptraj::Cluster::Control::AlgorithmArgs =
  "{hieragglo|dbscan|kmeans|dpeaks";

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

/** Set up clustering for a COORDS DataSet.*/
int Cpptraj::Cluster::Control::SetupForCoordsDataSet(DataSet_Coords* ds,
                                                     std::string const& maskExpr,
                                                     ArgList& analyzeArgs,
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
  mprintf("DEBUG: metric memory: %x\n", metric_);

  return Common(analyzeArgs);
}

/** Common setup. */
int Cpptraj::Cluster::Control::Common(ArgList& analyzeArgs) {
  // Allocate PairwiseMatrix.
  if (AllocatePairwise( analyzeArgs, metric_ )) {
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
  // Choose sieve restore option based on algorithm.
  if (sieve_ != 1) {
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
    //else if (bestRepStr == "cumulative_nosieve")
    //  bestRep_ = CUMULATIVE_NOSIEVE;
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


  Info();
  return 0;
}

// -----------------------------------------------------------------------------
void Cpptraj::Cluster::Control::Info() const {
  if (metric_    != 0) metric_->Info();
  if (algorithm_ != 0) algorithm_->Info();
}

int Cpptraj::Cluster::Control::Run() {
  if (metric_ == 0 || pmatrix_ == 0 || algorithm_ == 0) {
    mprinterr("Internal Error: Cluster::Control is not set up.\n");
    return 1;
  }
  // Set up the Metric
  if (metric_->Setup()) {
    mprinterr("Error: Metric setup failed.\n");
    return 1;
  }

  // Figure out which frames to cluster
  Sieve frameSieve;
  frameSieve.SetFramesToCluster(sieve_, metric_->Ntotal(), sieveSeed_);
  Cframes const& framesToCluster = frameSieve.FramesToCluster();

  // Cache distances if necessary
  pmatrix_->CacheDistances( framesToCluster );
  if (verbose_ > 1) pmatrix_->PrintCached();

  // Cluster
  if (algorithm_->DoClustering(clusters_, framesToCluster, *pmatrix_) != 0) {
    mprinterr("Error: Clustering failed.\n");
    return 1;
  }

  // ---------------------------------------------
  // Update cluster centroids here in case they need to be used to 
  // restore sieved frames
  clusters_.UpdateCentroids( metric_ );

  // Add sieved frames to existing clusters.
  if ( sieveRestore_ != NO_RESTORE ) {
    // Restore sieved frames
    mprintf("\tRestoring sieved frames.\n");
    switch (sieveRestore_) {
      case CLOSEST_CENTROID :
        clusters_.AddFramesByCentroid( frameSieve.SievedOut(), metric_ ); break;
      case EPSILON_CENTROID :
      case EPSILON_FRAME    :
        clusters_.AddFramesByCentroid( frameSieve.SievedOut(), metric_,
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

  // Find best representative frames for each cluster. TODO option to ignore sieved?
  if (BestReps::FindBestRepFrames(bestRep_, nRepsToSave_, clusters_, *pmatrix_, verbose_))
  {
    mprinterr("Error: Finding best representative frames for clusters failed.\n");
    return 1;
  }

  // DEBUG - print clusters to stdout
  if (verbose_ > 0) {
    mprintf("\nFINAL CLUSTERS:\n");
    clusters_.PrintClusters();
  }

  // TODO assign reference names

  

  return 0;
}

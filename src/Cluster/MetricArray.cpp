#include "MetricArray.h"
#include <algorithm> // std::min,max
#include <cmath> //sqrt
// For filling the pairwise cache
#ifdef _OPENMP
#include <omp.h>
#endif
#include "CentroidArray.h"
#include "Metric.h"
#include "../ArgList.h"
#include "../CpptrajFile.h"
#include "../CpptrajStdio.h"
#include "../DataFile.h"
#include "../DataFileList.h"
#include "../DataSet_PairwiseCache.h"
#include "../DataSetList.h"
#include "../OnlineVarT.h" // for metric stats average
#include "../ProgressBar.h"
#include "../StringRoutines.h"
// Metric classes
#include "Metric_RMS.h"
#include "Metric_DME.h"
#include "Metric_Scalar.h"
#include "Metric_SRMSD.h"
#include "Metric_Torsion.h"

/** CONSTRUCTOR */
Cpptraj::Cluster::MetricArray::MetricArray() :
  debug_(0),
  type_(MANHATTAN),
  ntotal_(0),
  cache_(0),
  cacheWasAllocated_(false),
  pw_mismatch_fatal_(true)
{}

/** DESTRUCTOR */
Cpptraj::Cluster::MetricArray::~MetricArray() {
  Clear();
}

/** COPY CONSTRUCTOR */
Cpptraj::Cluster::MetricArray::MetricArray(MetricArray const& rhs) :
  sets_(rhs.sets_),
  weights_(rhs.weights_),
  temp_(rhs.temp_),
  debug_(rhs.debug_),
  type_(rhs.type_),
  ntotal_(rhs.ntotal_),
  cache_(rhs.cache_),
  cacheWasAllocated_(rhs.cacheWasAllocated_),
  pw_mismatch_fatal_(rhs.pw_mismatch_fatal_)
{
  metrics_.reserve( rhs.metrics_.size() );
  for (std::vector<Metric*>::const_iterator it = rhs.metrics_.begin(); it != rhs.metrics_.end(); ++it)
    metrics_.push_back( (*it)->Copy() );
}

/** ASSIGNMENT */
Cpptraj::Cluster::MetricArray&
  Cpptraj::Cluster::MetricArray::operator=(MetricArray const& rhs)
{
  if (this == &rhs) return *this;
  metrics_.clear();
  metrics_.reserve( rhs.metrics_.size() );
  for (std::vector<Metric*>::const_iterator it = rhs.metrics_.begin(); it != rhs.metrics_.end(); ++it)
    metrics_.push_back( (*it)->Copy() );
  sets_ = rhs.sets_;
  weights_ = rhs.weights_;
  temp_ = rhs.temp_;
  debug_ = rhs.debug_;
  type_ = rhs.type_;
  ntotal_ = rhs.ntotal_;
  cache_ = rhs.cache_;
  cacheWasAllocated_ = rhs.cacheWasAllocated_;
  pw_mismatch_fatal_ = rhs.pw_mismatch_fatal_;
  return *this;
}

/** Clear the metric array. */
void Cpptraj::Cluster::MetricArray::Clear() {
  for (std::vector<Metric*>::iterator it = metrics_.begin(); it != metrics_.end(); ++it)
    delete *it;
  sets_.clear();
  weights_.clear();
  temp_.clear();
  cache_ = 0;
  cacheWasAllocated_ = false;
  pw_mismatch_fatal_ = true;
}

/** Pairwise args 1 */
const char* Cpptraj::Cluster::MetricArray::PairwiseArgs1_ =
  "[pairdist <name>] [pwrecalc]";

/** Pairwise args 2 */
const char* Cpptraj::Cluster::MetricArray::PairwiseArgs2_ =
  "[loadpairdist] [savepairdist] [pairwisecache {mem|disk|none}]";


/** Set up pairwise cache from arguments. */
int Cpptraj::Cluster::MetricArray::setupPairwiseCache(ArgList& analyzeArgs,
                                                      DataSetList& DSL,
                                                      DataFileList& DFL)
{
  if (metrics_.empty()) {
    mprinterr("Internal Error: allocatePairwise(): No Metrics.\n");
    return 1;
  }
  cacheWasAllocated_ = false;
  // The default pairwise cache file/dataset name
  const char* DEFAULT_PAIRDIST_NAME_ = "CpptrajPairDist";
  // The default pairwise cache file type
  DataFile::DataFormatType DEFAULT_PAIRDIST_TYPE_ =
#   ifdef BINTRAJ
    DataFile::CMATRIX_NETCDF;
#   else
    DataFile::CMATRIX_BINARY;
#   endif

  // Determine if we are saving/loading pairwise distances
  std::string pairdistname = analyzeArgs.GetStringKey("pairdist");
  DataFile::DataFormatType pairdisttype = DataFile::UNKNOWN_DATA;
  bool load_pair = analyzeArgs.hasKey("loadpairdist");
  bool save_pair = analyzeArgs.hasKey("savepairdist");
  // Check if we need to set a default file name
/*  std::string fname;
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
  }*/

  cache_ = 0;
  if (load_pair ||
      (!save_pair && !pairdistname.empty()))
  {
    // If 'loadpairdist' specified or 'pairdist' specified and 'savepairdist'
    // not specified, we either want to load from file or use an existing
    // data set.
    if (pairdistname.empty()) {
      pairdistname = DEFAULT_PAIRDIST_NAME_;
      pairdisttype = DEFAULT_PAIRDIST_TYPE_;
    }
    // First check if pairwise data exists
    DataSetList selected = DSL.SelectGroupSets( pairdistname, DataSet::PWCACHE );
    if (!selected.empty()) {
      if (selected.size() > 1)
        mprintf("Warning: '%s' matches multiple sets; only using '%s'\n",
                pairdistname.c_str(), selected[0]->legend());
      cache_ = (DataSet_PairwiseCache*)selected[0];
      mprintf("\tUsing existing pairwise set '%s'\n", cache_->legend());
    // Next check if file exists
    } else if (File::Exists( pairdistname )) {
      mprintf("\tLoading pairwise distances from file '%s'\n", pairdistname.c_str());
      DataFile dfIn;
      // TODO set data set name with ArgList?
      if (dfIn.ReadDataIn( pairdistname, ArgList(), DSL )) return 1;
      DataSet* ds = DSL.GetDataSet( pairdistname );
      if (ds == 0) return 1;
      if (ds->Group() != DataSet::PWCACHE) {
        mprinterr("Internal Error: AllocatePairwise(): Set is not a pairwise cache.\n");
        return 1;
      }
      cache_ = (DataSet_PairwiseCache*)ds;
      if (cache_ != 0 && save_pair) {
        mprintf("Warning: 'savepairdist' specified but pairwise cache loaded from file.\n"
                "Warning: Disabling 'savepairdist'.\n");
        save_pair = false;
      }
    } else
      pairdisttype = DEFAULT_PAIRDIST_TYPE_;

    if (cache_ == 0) {
      // Just 'pairdist' specified or loadpairdist specified and set/file not found.
      // Warn the user.
      mprintf("Warning: Pairwise distance matrix specified but cache/file '%s' not found.\n", pairdistname.c_str());
      if (!save_pair) {
        // If the file (or dataset) does not yet exist we will assume we want to save.
        mprintf("Warning: Pairwise distance matrix specified but not found; will save distances.\n");
        save_pair = true;
      }
    }
  } // END if load_pair

  // Create a pairwise cache if necessary
  if (cache_ == 0) {
    // Process DataSet type arguments
    DataSet::DataType pw_type = DataSet::PMATRIX_MEM;
    std::string pw_typeString = analyzeArgs.GetStringKey("pairwisecache");
    if (!pw_typeString.empty()) {
      if (pw_typeString == "mem")
        pw_type = DataSet::PMATRIX_MEM; 
      else if (pw_typeString == "disk") {
#       ifdef BINTRAJ
        pw_type = DataSet::PMATRIX_NC;
#       else
        // TODO regular disk file option
        mprinterr("Error: Pairwise disk cache requires NetCDF.\n");
        return 1;
#       endif
      } else if (pw_typeString == "none")
        pw_type = DataSet::UNKNOWN_DATA;
      else {
        mprinterr("Error: Unrecognized option for 'pairwisecache' ('%s')\n", pw_typeString.c_str());
        return 1;
      }
    }
    // Allocate cache if necessary
    if (pw_type != DataSet::UNKNOWN_DATA) {
      MetaData meta;
      if (!pairdistname.empty())
        meta.SetName( pairdistname );
      else
        meta.SetName( DSL.GenerateDefaultName("CMATRIX") );
      // Cache-specific setup.
      if (pw_type == DataSet::PMATRIX_NC)
        meta.SetFileName( pairdistname ); // TODO separate file name?
      //cache_ = (DataSet_PairwiseCache*)DSL.AddSet( pw_type, meta, "CMATRIX" );
      // To maintain compatibility with pytraj, the cluster number vs time
      // set **MUST** be allocated before the cache. Set up outside the
      // DataSetList here and set up; add to DataSetList later in Setup.
      cache_ = (DataSet_PairwiseCache*)DSL.AllocateSet( pw_type, meta );
      if (cache_ == 0) {
        mprinterr("Error: Could not allocate pairwise cache.\n");
        return 1;
      }
      cacheWasAllocated_ = true;
      if (debug_ > 0)
        mprintf("DEBUG: Allocated pairwise distance cache: %s\n", cache_->legend());
    }
  }

  // Setup pairwise matrix
  //if (pmatrix_.Setup(metric_, cache_)) return 1;

  if (save_pair) {
    if (cache_ == 0) {
      mprintf("Warning: Not caching distances; ignoring 'savepairdist'\n");
    } else {
      if (pairdistname.empty())
        pairdistname = DEFAULT_PAIRDIST_NAME_;
      if (pairdisttype == DataFile::UNKNOWN_DATA)
        pairdisttype = DEFAULT_PAIRDIST_TYPE_;
      // TODO enable saving for other set types?
      if (cache_->Type() == DataSet::PMATRIX_MEM) {
        DataFile* pwd_file = DFL.AddDataFile( pairdistname, pairdisttype, ArgList() );
        if (pwd_file == 0) return 1;
        pwd_file->AddDataSet( cache_ );
        if (debug_ > 0)
          mprintf("DEBUG: Saving pw distance cache '%s' to file '%s'\n", cache_->legend(),
                  pwd_file->DataFilename().full());
      }
    }
  }
  pw_mismatch_fatal_ = !analyzeArgs.hasKey("pwrecalc");

  return 0;
}

/** \return Pointer to first Metric related to COORDS, or 0 if not present. */
Cpptraj::Cluster::Metric const*
  Cpptraj::Cluster::MetricArray::CoordsMetric() const
{
  for (std::vector<Metric*>::const_iterator it = metrics_.begin(); it != metrics_.end(); ++it)
  {
    if ( (*it)->MetricType() == Metric::RMS ||
         (*it)->MetricType() == Metric::DME ||
         (*it)->MetricType() == Metric::SRMSD )
      return *it;
  }
  return 0;
}

/** /return Pointer to Metric of given type. */
Cpptraj::Cluster::Metric* Cpptraj::Cluster::MetricArray::AllocateMetric(Metric::Type mtype)
{
  Metric* met = 0;
  switch (mtype) {
    case Metric::RMS       : met = new Metric_RMS(); break;
    case Metric::DME       : met = new Metric_DME(); break;
    case Metric::SRMSD     : met = new Metric_SRMSD(); break;
    case Metric::SCALAR    : met = new Metric_Scalar(); break;
    case Metric::TORSION   : met = new Metric_Torsion(); break;
    default: mprinterr("Error: Unhandled Metric in AllocateMetric.\n");
  }
  return met;
}

/** Recognized metric args. */
const char* Cpptraj::Cluster::MetricArray::MetricArgs_ =
  "[{dme|rms|srmsd} [mass] [nofit] [<mask>]] [{euclid|manhattan}] [wgt <list>]";

/** Initialize with given sets and arguments. */
int Cpptraj::Cluster::MetricArray::initMetricArray(DataSetList const& setsToCluster, ArgList& analyzeArgs)
{
  // Get rid of any previous metrics
  Clear();
  // Get arguments for any COORDS metrics
  int usedme = (int)analyzeArgs.hasKey("dme");
  int userms = (int)analyzeArgs.hasKey("rms");
  int usesrms = (int)analyzeArgs.hasKey("srmsd");
  bool useMass = analyzeArgs.hasKey("mass");
  bool nofit   = analyzeArgs.hasKey("nofit");
  std::string maskExpr = analyzeArgs.GetMaskNext();
  // Get other args
  if (analyzeArgs.hasKey("euclid"))
    type_ = EUCLID;
  else if (analyzeArgs.hasKey("manhattan"))
    type_ = MANHATTAN;
  else {
    // Default
    if (setsToCluster.size() > 1)
      type_ = EUCLID;
    else
      type_ = MANHATTAN;
  }
  std::string wgtArgStr = analyzeArgs.GetStringKey("wgt");

  // Check args
  if (usedme + userms + usesrms > 1) {
    mprinterr("Error: Specify either 'dme', 'rms', or 'srmsd'.\n");
    return 1;
  }
  Metric::Type coordsMetricType;
  if      (usedme)  coordsMetricType = Metric::DME;
  else if (userms)  coordsMetricType = Metric::RMS;
  else if (usesrms) coordsMetricType = Metric::SRMSD;
  else coordsMetricType = Metric::RMS; // default

  // For each input set, set up the appropriate metric
  for (DataSetList::const_iterator ds = setsToCluster.begin(); ds != setsToCluster.end(); ++ds)
  {
    Metric::Type mtype = Metric::UNKNOWN_METRIC;
    if ( (*ds)->Group() == DataSet::COORDINATES )
    {
      mtype = coordsMetricType;
    } else if ((*ds)->Group() == DataSet::SCALAR_1D) {
      if ( (*ds)->Meta().IsTorsionArray() )
        mtype = Metric::TORSION;
      else
        mtype = Metric::SCALAR;
    }

    if (mtype == Metric::UNKNOWN_METRIC) {
      mprinterr("Error: Set '%s' is a type not yet supported by Cluster.\n", (*ds)->legend());
      return 1;
    }

    Metric* met = AllocateMetric( mtype );
    if (met == 0) {
      mprinterr("Internal Error: Could not allocate metric for set '%s'\n", (*ds)->legend());
      return 1;
    }

    // Initialize the metric
    int err = 0;
    switch (mtype) {
      case Metric::RMS :
        err = ((Metric_RMS*)met)->Init((DataSet_Coords*)*ds, AtomMask(maskExpr), nofit, useMass); break;
      case Metric::DME :
        err = ((Metric_DME*)met)->Init((DataSet_Coords*)*ds, AtomMask(maskExpr)); break;
      case Metric::SRMSD :
        err = ((Metric_SRMSD*)met)->Init((DataSet_Coords*)*ds, AtomMask(maskExpr), nofit, useMass, debug_); break;
      case Metric::SCALAR :
        err = ((Metric_Scalar*)met)->Init((DataSet_1D*)*ds); break;
      case Metric::TORSION :
        err = ((Metric_Torsion*)met)->Init((DataSet_1D*)*ds); break;
      default:
        mprinterr("Error: Unhandled Metric setup.\n");
        err = 1;
    }
    if (err != 0) {
      mprinterr("Error: Metric setup failed.\n");
      return 1;
    }
    metrics_.push_back( met );
    sets_.push_back( *ds );

  } // END loop over input data sets

  // Process weight args if needed
  if (!wgtArgStr.empty()) {
    ArgList wgtArgs(wgtArgStr, ",");
    // Need 1 arg for every set 
    if (wgtArgs.Nargs() != (int)metrics_.size()) {
      mprinterr("Error: Expected %zu comma-separated args for wgt, got %i\n",
                metrics_.size(), wgtArgs.Nargs());
      return 1;
    }
    weights_.reserve( wgtArgs.Nargs() );
    for (int arg = 0; arg != wgtArgs.Nargs(); arg++)
      weights_.push_back( convertToDouble( wgtArgs[arg] ) );
  } else {
    // Default weights are 1
    weights_.assign(metrics_.size(), 1.0);
  }
  // Temp space for calcs
  temp_.assign(metrics_.size(), 0.0);

  return 0;
}

/** Initialize metrics and pairwise cache (if needed). */
int Cpptraj::Cluster::MetricArray::Initialize(DataSetList const& setsToCluster,
                                              DataSetList& dslIn, DataFileList& dflIn,
                                              ArgList& analyzeArgs, int debugIn)
{
  debug_ = debugIn;
  if (initMetricArray(setsToCluster, analyzeArgs)) {
    mprinterr("Error: Metric initialization failed.\n");
    return 1;
  }
  if (setupPairwiseCache(analyzeArgs, dslIn, dflIn)) {
    mprinterr("Error: Pairwise cache initialization failed.\n");
    return 1;
  }
  return 0;
}

/** Call the Setup function for all metrics. Check that size of each Metric
  * is the same.
  */
int Cpptraj::Cluster::MetricArray::Setup() {
  int err = 0;
  ntotal_ = 0;
  for (std::vector<Metric*>::const_iterator it = metrics_.begin();
                                            it != metrics_.end(); ++it)
  {
    if ((*it)->Setup()) {
      mprinterr("Error: Metric '%s' setup failed.\n", (*it)->Description().c_str());
      err++;
    }
    if (ntotal_ == 0)
      ntotal_ = (*it)->Ntotal();
    else {
      if ( (*it)->Ntotal() != ntotal_ ) {
        mprinterr("Error: Number of points covered by metric '%s' (%u) is not equal\n"
                  "Error:  to number of points covered by metric '%s' (%u)\n",
                  (*it)->Description().c_str(), (*it)->Ntotal(),
                  metrics_.front()->Description().c_str(), metrics_.front()->Ntotal());
        return 1;
      }
    }
  }
  return err;
}

/** Call the Info function for all metrics. */
void Cpptraj::Cluster::MetricArray::Info() const {
  static const char* DistanceTypeStr[] = { "Manhattan", "Euclidean" };
  mprintf("\tUsing %s distance.\n", DistanceTypeStr[type_]);
  for (unsigned int idx = 0; idx != metrics_.size(); idx++)
  {
    mprintf("\tMetric %u for '%s', weight factor %g\n",
            idx, sets_[idx]->legend(), weights_[idx]);
    metrics_[idx]->Info();
  }

  if (cache_ == 0)
    mprintf("\tPairwise distances will not be cached.\n");
  else {
    if (cache_->Size() > 0)
      mprintf("\tUsing existing pairwise cache: %s (%s)\n",
              cache_->legend(), cache_->description());
    else
      mprintf("\tPairwise distances will be cached: %s (%s)\n",
              cache_->legend(), cache_->description());
    if (pw_mismatch_fatal_)
      mprintf("\tCalculation will be halted if frames in cache do not match.\n");
    else
      mprintf("\tPairwise distances will be recalculated if frames in cache do not match.\n");
  }
}

// -----------------------------------------------
/** Calculate new centroids for given list. */
void Cpptraj::Cluster::MetricArray::NewCentroid(CentroidArray& centroids, Cframes const& framesIn)
const
{
  //mprintf("DEBUG: Calling MetricArray::NewCentroid for %zu frames.\n", framesIn.size());
  centroids.Clear();
  for (std::vector<Metric*>::const_iterator it = metrics_.begin(); it != metrics_.end(); ++it)
  {
    centroids.push_back( (*it)->NewCentroid( framesIn ) );
  }
  //mprintf("DEBUG: There are %zu centroids for %zu metrics.\n", centroids.size(), metrics_.size());
}

/** Calculate centroids in given list. */
void Cpptraj::Cluster::MetricArray::CalculateCentroid(CentroidArray& centroids, Cframes const& framesIn)
const
{
  if (centroids.size() != metrics_.size()) {
    mprinterr("Internal Error: MetricArray::CalculateCentroid: centroids (%u) and metrics_ (%zu) sizes do not match.\n", centroids.size(), metrics_.size());
    return;
  }
  for (unsigned int idx = 0; idx != metrics_.size(); idx++)
    metrics_[idx]->CalculateCentroid( centroids[idx], framesIn );
}

/** Update centroids by performing given operation between given frame and centroids. */
void Cpptraj::Cluster::MetricArray::FrameOpCentroid(int f1, CentroidArray& centroids,
                                                    double oldSize, Metric::CentOpType OP)
const
{
  if (centroids.size() != metrics_.size()) {
    mprinterr("Internal Error: MetricArray::FrameOpCentroid: centroids (%u) and metrics_ (%zu) sizes do not match.\n", centroids.size(), metrics_.size());
    return;
  }
  for (unsigned int idx = 0; idx != metrics_.size(); idx++)
    metrics_[idx]->FrameOpCentroid( f1, centroids[idx], oldSize, OP );
}

/** Manhattan distance. */
double Cpptraj::Cluster::MetricArray::Dist_Manhattan(std::vector<double> const& arrayIn) const {
  double dist = 0.0;
  for (unsigned int idx = 0; idx != arrayIn.size(); idx++)
    dist += (weights_[idx] * arrayIn[idx]);
  return dist;
}

/** Euclidean distance. */
double Cpptraj::Cluster::MetricArray::Dist_Euclidean(std::vector<double> const& arrayIn) const {
  double sumdist2 = 0.0;
  for (unsigned int idx = 0; idx != arrayIn.size(); idx++) {
    double dist2 = arrayIn[idx] * arrayIn[idx];
    sumdist2 += (weights_[idx] * dist2);
  }
  return sqrt(sumdist2);
} 

/** Perform distance calc according to current type. */
double Cpptraj::Cluster::MetricArray::DistCalc(std::vector<double> const& arrayIn) const {
  double dist = 0;
  switch (type_) {
    case MANHATTAN : dist = Dist_Manhattan(arrayIn); break;
    case EUCLID    : dist = Dist_Euclidean(arrayIn); break;
  }
  return dist;
}

/** Calculate distance between given frame and centroids. */
double Cpptraj::Cluster::MetricArray::FrameCentroidDist(int frame, CentroidArray const& centroids )
{
  if (centroids.size() != metrics_.size()) {
    mprinterr("Internal Error: MetricArray::FrameCentroidDist: centroids (%u) and metrics_ (%zu) sizes do not match.\n", centroids.size(), metrics_.size());
    return -1.0;
  }
  for (unsigned int idx = 0; idx != metrics_.size(); idx++)
    temp_[idx] = metrics_[idx]->FrameCentroidDist( frame, centroids[idx] );
  return DistCalc(temp_);
}

/** Calculate distance between given centroids. */
double Cpptraj::Cluster::MetricArray::CentroidDist(CentroidArray const& C1, CentroidArray const& C2)
{
  if (C1.size() != metrics_.size() || C2.size() != metrics_.size()) {
    mprinterr("Internal Error: MetricArray::CentroidDist: centroids (%u, %u) and metrics_ (%zu) sizes do not match.\n", C1.size(), C2.size(), metrics_.size());
    return -1.0;
  }
  for (unsigned int idx = 0; idx != metrics_.size(); idx++)
    temp_[idx] = metrics_[idx]->CentroidDist( C1[idx], C2[idx] );
  return DistCalc(temp_);
}

/** \return distance between frames (uncached). */
double Cpptraj::Cluster::MetricArray::Uncached_Frame_Distance(int f1, int f2)
{
  for (unsigned int idx = 0; idx != metrics_.size(); idx++)
    temp_[idx] = metrics_[idx]->FrameDist(f1, f2);
  return DistCalc(temp_);
}

/** \return distance between frames (cached or uncached). */
double Cpptraj::Cluster::MetricArray::Frame_Distance(int f1, int f2) {
  if (cache_ != 0)
  {
    // TODO protect against f1/f2 out of bounds.
    int idx1 = cache_->FrameToIdx()[f1];
    if (idx1 != -1) {
      int idx2 = cache_->FrameToIdx()[f2];
      if (idx2 != -1)
        return cache_->CachedDistance(idx1, idx2);
    }
  }
  // If here, distance was not cached or no cache.
  return Uncached_Frame_Distance(f1, f2); 
}

// -----------------------------------------------
/** Loop over all pairs. For each distance, evaluate the contribution
  * of each metric to the total distance (Manhattan).
  */
void Cpptraj::Cluster::MetricArray::CalculateMetricContributions(Cframes const& framesToCluster,
                                                                 CpptrajFile& outfile)
{
  if (metrics_.size() < 2) {
    mprintf("\tLess than 2 metrics; skipping metric contribution calculation.\n");
    return;
  }
  mprintf("\tCalculating metric contributions to total distance:\n");
  // This array will accumulate contribution fractions
  std::vector<Stats<double>> mfrac( metrics_.size() );
  // This array will accumulate contribution averages
  std::vector<Stats<double>> mavg( metrics_.size() );
  // This array will hold each contribution, adjust for distance type and weight
  std::vector<double> mcont( metrics_.size() );
  // These arrays will record the minimum and maximum distance contributions
  std::vector<double> mmin;
  std::vector<double> mmax;
//  mprintf("DEBUG: %8s %8s %8s %8s", "Frame1", "Frame2", "Distance", "Sum"); // DEBUG
//  for (unsigned int idx = 0; idx != metrics_.size(); idx++) // DEBUG
//    mprintf(" {%8s %8s %6s}", "Raw", "Adj.", "Frac"); // DEBUG
//  mprintf("\n"); // DEBUG
  for (Cframes_it frm1 = framesToCluster.begin(); frm1 != framesToCluster.end(); ++frm1)
  {
    for (Cframes_it frm2 = frm1 + 1; frm2 != framesToCluster.end(); ++frm2)
    {
      // Populate the temp array
      //double dist = Uncached_Frame_Distance(*frm1, *frm2); // DEBUG
      Uncached_Frame_Distance(*frm1, *frm2);
      // Do individual contributions, adjusted for weights and distance type
      double sum = 0;
      if (type_ == MANHATTAN) {
        for (unsigned int idx = 0; idx != temp_.size(); idx++) {
          mcont[idx] = (weights_[idx] * temp_[idx]);
          sum += mcont[idx];
        }
      } else if (type_ == EUCLID) {
        for (unsigned int idx = 0; idx != temp_.size(); idx++) {
          mcont[idx] = (weights_[idx] * (temp_[idx]*temp_[idx]));
          sum += mcont[idx];
        }
      } else {
        // Sanity check
        mprinterr("Internal Error: CalculateMetricContributions: Unhandled metric summation.\n");
        return;
      }
      // Do fraction contributions
      for (unsigned int idx = 0; idx != temp_.size(); idx++)
        mfrac[idx].accumulate( mcont[idx] / sum );
      // Set initial min/max values
      if (mmin.empty()) {
        mmin.resize( metrics_.size() );
        mmax.resize( metrics_.size() );
        for (unsigned int idx = 0; idx != temp_.size(); idx++) {
          mmin[idx] = temp_[idx];
          mmax[idx] = temp_[idx];
        }
      }
      // Do averages
      //mprintf("DEBUG: %8i %8i %8.3f %8.3f", *frm1 + 1, *frm2 + 1, dist, sum); // DEBUG
      for (unsigned int idx = 0; idx != temp_.size(); idx++)
      {
        mavg[idx].accumulate( temp_[idx] );
        mmin[idx] = std::min( temp_[idx], mmin[idx] );
        mmax[idx] = std::max( temp_[idx], mmax[idx] );
        //mprintf(" {%8.3f %8.3f %6.2f}", temp_[idx], mcont[idx], mcont[idx] / sum); // DEBUG
      }
      //mprintf("\n"); // DEBUG
    }
  }
  // Output
  outfile.Printf("%-7s %6s %6s %12s %12s %12s %12s %s\n",
                 "#Metric", "FracAv", "FracSD", "Avg", "SD", "Min", "Max", "Description");
  for (unsigned int idx = 0; idx != mfrac.size(); idx++) {
    outfile.Printf("%-7u %6.4f %6.4f %12.4f %12.4f %12.4f %12.4f \"%s\"\n", idx,
                   mfrac[idx].mean(), sqrt(mfrac[idx].variance()),
                   mavg[idx].mean(), sqrt(mavg[idx].variance()),
                   mmin[idx], mmax[idx],
                   metrics_[idx]->Description().c_str());
  }
  outfile.Flush();
}

// -----------------------------------------------
/** \return A string containing the description of all metrics. */
std::string Cpptraj::Cluster::MetricArray::descriptionOfMetrics() const {
  std::string out;
  for (std::vector<Metric*>::const_iterator it = metrics_.begin();
                                            it != metrics_.end(); ++it)
  {
    if (it != metrics_.begin())
      out.append(",");
    out.append( (*it)->Description() );
  }
  return out;
}

/** Request that distances for frames in given array be cached.
  * \param framesToCache the frames to cache.
  * \param sieveIn Sieve value (if any) used to generate framesToCache. This is
  *        purely for bookkeeping inside DataSet_PairwiseCache.
  * When pw_mismatch_fatal_ is true, if a cache is present but the frames to
  * cache do not match what is in the cache, exit with an error;
  * otherwise recalculate the cache.
  */
int Cpptraj::Cluster::MetricArray::CacheDistances(Cframes const& framesToCache, int sieveIn)
{
  if (framesToCache.size() < 1) return 0;
  // If no cache we can leave.
  if (cache_ == 0) return 0;

  bool do_cache = true;
  if (cache_->Size() > 0) {
    do_cache = false;
    mprintf("\tUsing existing cache '%s'\n", cache_->legend());
    // If cache is already populated, check that it is valid.
    // The frames to cache must match cached frames.
    if (!cache_->CachedFramesMatch( framesToCache )) {
      if (pw_mismatch_fatal_) {
        mprinterr("Error: Frames to cache do not match those in existing cache.\n");
        return 1;
      } else {
        mprintf("Warning: Frames to cache do not match those in existing cache.\n"
                "Warning: Re-calculating pairwise cache.\n");
        do_cache = true;
      }
    }
    // TODO Check metric? Total frames?
  }

  if (do_cache) {
    // Sanity check
    if (metrics_.empty()) {
      mprinterr("Internal Error: MetricArray::CacheDistances(): Metric is null.\n");
      return 1;
    }
    // Cache setup
    if (cache_->SetupCache( Ntotal(), framesToCache, sieveIn, descriptionOfMetrics() ))
      return 1;
    // Fill cache
    if (calcFrameDistances(framesToCache))
      return 1;
  }
  return 0;
}

/** Cache distances between given frames using SetElement(). */
int Cpptraj::Cluster::MetricArray::calcFrameDistances(Cframes const& framesToCache)
{
  mprintf("\tCaching distances for %zu frames.\n", framesToCache.size());

  int f2end = (int)framesToCache.size();
  int f1end = f2end - 1;
  ParallelProgress progress(f1end);
  int f1, f2;
  // For OMP, every other thread will need its own Cdist.
  MetricArray* MyMetrics = this;
# ifdef _OPENMP
# pragma omp parallel private(MyMetrics, f1, f2) firstprivate(progress)
  {
  int mythread = omp_get_thread_num();
  progress.SetThread( mythread );
  if (mythread == 0) {
    mprintf("\tParallelizing pairwise distance calc with %i threads\n", omp_get_num_threads());
    MyMetrics = this;
  } else
    MyMetrics = new MetricArray(*this);
# pragma omp for schedule(dynamic)
# endif
  for (f1 = 0; f1 < f1end; f1++) {
    progress.Update(f1);
    for (f2 = f1 + 1; f2 < f2end; f2++) {
      cache_->SetElement( f1, f2, MyMetrics->Uncached_Frame_Distance(framesToCache[f1], framesToCache[f2]) );
    }
  }
# ifdef _OPENMP
  if (mythread > 0)
    delete MyMetrics;
  } // END omp parallel
# endif
  progress.Finish();
  return 0;
}

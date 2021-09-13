#include "MetricArray.h"
#include <cmath> //sqrt
#include "CentroidArray.h"
#include "Metric.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../DataSetList.h"
#include "../StringRoutines.h"
// Metric classes
#include "Metric_RMS.h"
#include "Metric_DME.h"
#include "Metric_Scalar.h"
#include "Metric_SRMSD.h"
#include "Metric_Torsion.h"

/** CONSTRUCTOR */
Cpptraj::Cluster::MetricArray::MetricArray() :
  type_(MANHATTAN),
  ntotal_(0)
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
  type_(rhs.type_),
  ntotal_(rhs.ntotal_)
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
  type_ = rhs.type_;
  ntotal_ = rhs.ntotal_;
  return *this;
}

/** Clear the metric array. */
void Cpptraj::Cluster::MetricArray::Clear() {
  for (std::vector<Metric*>::iterator it = metrics_.begin(); it != metrics_.end(); ++it)
    delete *it;
  sets_.clear();
  weights_.clear();
  temp_.clear();
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

/** Recognized args. */
const char* Cpptraj::Cluster::MetricArray::MetricArgs_ =
  "[{dme|rms|srmsd} [mass] [nofit] [<mask>]] [{euclid|manhattan}] [wgt <list>]";

/** Initialize with given sets and arguments. */
int Cpptraj::Cluster::MetricArray::InitMetricArray(DataSetList const& dslIn, ArgList& analyzeArgs,
                                                   int debugIn)
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
    if (dslIn.size() > 1)
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
  for (DataSetList::const_iterator ds = dslIn.begin(); ds != dslIn.end(); ++ds)
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
        err = ((Metric_SRMSD*)met)->Init((DataSet_Coords*)*ds, AtomMask(maskExpr), nofit, useMass, debugIn); break;
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
}

// -----------------------------------------------
/** Calculate new centroids for given list. */
void Cpptraj::Cluster::MetricArray::NewCentroid(CentroidArray& centroids, Cframes const& framesIn)
// TODO const?
{
  centroids.Clear();
  for (std::vector<Metric*>::iterator it = metrics_.begin(); it != metrics_.end(); ++it)
  {
    centroids.push_back( (*it)->NewCentroid( framesIn ) );
  }
}

/** Calculate centroids in given list. */
void Cpptraj::Cluster::MetricArray::CalculateCentroid(CentroidArray& centroids, Cframes const& framesIn)
{
  if (centroids.size() != metrics_.size()) {
    mprinterr("Internal Error: MetricArray::CalculateCentroid: centroids and metrics_ sizes do not match.\n");
    return;
  }
  for (unsigned int idx = 0; idx != metrics_.size(); idx++)
    metrics_[idx]->CalculateCentroid( centroids[idx], framesIn );
}

/** Update centroids by performing given operation between given frame and centroids. */
void Cpptraj::Cluster::MetricArray::FrameOpCentroid(int f1, CentroidArray& centroids,
                                                    double oldSize, Metric::CentOpType OP)
{
  if (centroids.size() != metrics_.size()) {
    mprinterr("Internal Error: MetricArray::FrameOpCentroid: centroids and metrics_ sizes do not match.\n");
    return;
  }
  for (unsigned int idx = 0; idx != metrics_.size(); idx++)
    metrics_[idx]->FrameOpCentroid( f1, centroids[idx], oldSize, OP );
}

/** Manhattan distance. */
double Cpptraj::Cluster::MetricArray::Dist_Manhattan() const {
  double dist = 0.0;
  for (unsigned int idx = 0; idx != temp_.size(); idx++)
    dist += (weights_[idx] * temp_[idx]);
  return dist;
}

/** Euclidean distance. */
double Cpptraj::Cluster::MetricArray::Dist_Euclidean() const {
  double sumdist2 = 0.0;
  for (unsigned int idx = 0; idx != temp_.size(); idx++) {
    double dist2 = temp_[idx] * temp_[idx];
    sumdist2 += (weights_[idx] * dist2);
  }
  return sqrt(sumdist2);
} 

/** Perform distance calc according to current type. */
double Cpptraj::Cluster::MetricArray::DistCalc() const {
  double dist = 0;
  switch (type_) {
    case MANHATTAN : dist = Dist_Manhattan(); break;
    case EUCLID    : dist = Dist_Euclidean(); break;
  }
  return dist;
}

/** Calculate distance between given frame and centroids. */
double Cpptraj::Cluster::MetricArray::FrameCentroidDist(int frame, CentroidArray const& centroids )
{
  if (centroids.size() != metrics_.size()) {
    mprinterr("Internal Error: MetricArray::FrameOpCentroid: centroids and metrics_ sizes do not match.\n");
    return -1.0;
  }
  for (unsigned int idx = 0; idx != metrics_.size(); idx++)
    temp_[idx] = metrics_[idx]->FrameCentroidDist( frame, centroids[idx] );
  return DistCalc();
}

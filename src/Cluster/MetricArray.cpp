#include "MetricArray.h"
#include "Metric.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../DataSetList.h"
// Metric classes
#include "Metric_RMS.h"
#include "Metric_DME.h"
#include "Metric_Scalar.h"
#include "Metric_SRMSD.h"
#include "Metric_Torsion.h"

/** CONSTRUCTOR */
Cpptraj::Cluster::MetricArray::MetricArray() {}

/** DESTRUCTOR */
Cpptraj::Cluster::MetricArray::~MetricArray() {
  Clear();
}

/** COPY CONSTRUCTOR */
Cpptraj::Cluster::MetricArray::MetricArray(MetricArray const& rhs) {
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
  return *this;
}

/** Clear the metric array. */
void Cpptraj::Cluster::MetricArray::Clear() {
  for (std::vector<Metric*>::iterator it = metrics_.begin(); it != metrics_.end(); ++it)
    delete *it;
  sets_.clear();
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

  } // END loop over input data sets

  return 0;
} 

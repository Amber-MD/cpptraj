#include "Control.h"
#include "../CpptrajStdio.h"
#include "../DataSet_Coords.h"
// PairwiseMatrix classes
#include "PairwiseMatrix_MEM.h"
// Metric classes
#include "Metric_RMS.h"
// Algorithms
#include "Algorithm_HierAgglo.h"

void Cpptraj::Cluster::Control::Help() {
  mprintf("[crdset <COORDS set>]\n");
}

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
Cpptraj::Cluster::Algorithm* Cpptraj::Cluster::Control::AllocateAlgorithm(Algorithm::Type atype)
{
  Algorithm* alg = 0;
  switch (atype) {
    case Algorithm::HIERAGGLO : alg = new Algorithm_HierAgglo(); break;
    default : mprinterr("Error: Unhandled Algorithm in AllocateAlgorithm.\n");
  }
  return alg;
}

/** Set up Algorithm from keyword + arguments. */
int Cpptraj::Cluster::Control::AllocateAlgorithm(ArgList& analyzeArgs) {
  Algorithm::Type atype;
  if      (analyzeArgs.hasKey("hieragglo")) atype = Algorithm::HIERAGGLO;
  else if (analyzeArgs.hasKey("dbscan"   )) atype = Algorithm::DBSCAN;
  else if (analyzeArgs.hasKey("kmeans") ||
           analyzeArgs.hasKey("means")    ) atype = Algorithm::KMEANS;
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
                                                     int verbose)
{
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
      err = ((Metric_RMS*)metric_)->Setup(ds, AtomMask(maskExpr), nofit, useMass); break;
    default:
      mprinterr("Error: Unhandled Metric setup.\n");
      err = 1;
  }
  if (err != 0) {
    mprinterr("Error: Metric setup failed.\n");
    return 1;
  }

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

  return 0;
}

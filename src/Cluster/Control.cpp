#include "Control.h"
#include "../CpptrajStdio.h"
#include "../DataSet_Coords.h"
// PairwiseMatrix classes
#include "PairwiseMatrix_MEM.h"

void Cpptraj::Cluster::Control::Help() {
  mprintf("[crdset <COORDS set>]\n");
}

int Cpptraj::Cluster::Control::allocatePairwise(ArgList& analyzeArgs, Metric::Type mtype)
{
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
  pmatrix_ = 0;
  switch (pw_type) {
    case PairwiseMatrix::MEM : pmatrix_ = new PairwiseMatrix_MEM(mtype); break;
    default: mprinterr("Error: Unhandled pairwise type.\n");
  }
  if (pmatrix_ == 0) return 1;
  return 0;
}


int Cpptraj::Cluster::Control::SetupForCoordsDataSet(DataSet_Coords* ds, ArgList& analyzeArgs,
                                                     int verbose)
{
  // Determine metric. Valid ones for COORDS are RMS, DMR, SRMSD
  
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
  if (allocatePairwise( analyzeArgs, mtype )) return 1;

  return 0;
}

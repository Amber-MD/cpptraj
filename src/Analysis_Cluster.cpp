#include "Analysis_Cluster.h"
#include "CpptrajStdio.h"

using namespace Cpptraj::Cluster;

// Analysis_Cluster::Help()
void Analysis_Cluster::Help() const {

}

// Analysis_Cluster::Setup()
Analysis::RetType Analysis_Cluster::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringKey("crdset");
  coords_ = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
  if (coords_ == 0) {
    mprinterr("Error: Could not locate COORDS set corresponding to %s\n",
              setname.c_str());
    return Analysis::ERR;
  }

  // Get the mask string 
  std::string maskexpr = analyzeArgs.GetMaskNext();
  /*if (!refs_.empty() && refmaskexpr_.empty()) {
    refmaskexpr_ = maskexpr_;
    if (refmaskexpr_.empty()) {
      refmaskexpr_.assign("!@H=");
      mprintf("Warning: 'assignrefs' specified but no 'refmask' given.\n"
              "Warning:   Using default mask expression: '%s'\n", refmaskexpr_.c_str());
    }
  }*/


  control_.SetupForCoordsDataSet(coords_, maskexpr, analyzeArgs, debugIn);

  return Analysis::OK;

}

// Analysis_Cluster::Analyze()
Analysis::RetType Analysis_Cluster::Analyze() {
  int err = control_.Run();
  if (err != 0)
    return Analysis::ERR;
  return Analysis::OK;
}

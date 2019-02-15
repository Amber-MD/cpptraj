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

  if (control_.SetupForCoordsDataSet(coords_, analyzeArgs, setup.DSL(), setup.DFL(), debugIn))
    return Analysis::ERR;

  control_.Info();

  masterDSL_ = setup.DslPtr();

  return Analysis::OK;

}

// Analysis_Cluster::Analyze()
Analysis::RetType Analysis_Cluster::Analyze() {
  Timer t_total;
  t_total.Start();
  if (control_.Run() != 0) {
    mprinterr("Error: Clustering failed.\n");
    return Analysis::ERR;
  }
  // Output
  if (control_.Output(*masterDSL_) != 0) {
    mprinterr("Error: Cluster output failed.\n");
    return Analysis::ERR;
  }
  t_total.Stop();
  control_.Timing(t_total.Total());
  return Analysis::OK;
}

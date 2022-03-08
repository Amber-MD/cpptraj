#include "Analysis_Clustering.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords.h"

using namespace Cpptraj::Cluster;

// Analysis_Clustering::Help()
void Analysis_Clustering::Help() const {
  mprintf("\t[crdset <crd set>] [data <dset0>[,<dset1>...]] [nocoords]\n");
  Control::Help();
}

// Analysis_Clustering::Setup()
Analysis::RetType Analysis_Clustering::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  control_.SetDebug( debugIn );
  DataSetList inputDsets;
  DataSet_Coords* coords = 0;
  // First check for data
  std::string dataSetname = analyzeArgs.GetStringKey("data");
  if (!dataSetname.empty()) {
    // Get data sets for clustering
    ArgList dsnames(dataSetname, ",");
    for (ArgList::const_iterator name = dsnames.begin(); name != dsnames.end(); ++name) {
      DataSetList tempDSL = setup.DSL().GetMultipleSets( *name );
      if (tempDSL.empty()) {
        mprinterr("Error: %s did not correspond to any data sets.\n", dataSetname.c_str());
        return Analysis::ERR;
      }
      inputDsets += tempDSL;
    }
  }
  // If 'crdset' specified or if 'nocoords' *not* specified, attempt to get
  // COORDS DataSet from master DataSetList.
  std::string setname = analyzeArgs.GetStringKey("crdset");
  if (!setname.empty() || !analyzeArgs.hasKey("nocoords")) {
    coords = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
    if (coords == 0) {
      mprinterr("Error: Could not locate COORDS set corresponding to %s\n",
                setname.c_str());
      return Analysis::ERR;
    }
    // If 'data' was not specified, this is the set we will cluster on.
    if (dataSetname.empty())
      inputDsets.AddCopyOfSet( coords );
  }

  if (control_.SetupClustering(inputDsets, coords, analyzeArgs, setup.DSL(), setup.DFL(), debugIn))
      return Analysis::ERR;

  control_.Info();

  masterDSL_ = setup.DslPtr();

  return Analysis::OK;

}

// Analysis_Clustering::Analyze()
Analysis::RetType Analysis_Clustering::Analyze() {
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

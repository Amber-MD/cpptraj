#include "Analysis_Cluster.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords.h"
#include "Cluster/Metric_Data.h"

using namespace Cpptraj::Cluster;

// Analysis_Cluster::Help()
void Analysis_Cluster::Help() const {

}

// Analysis_Cluster::Setup()
Analysis::RetType Analysis_Cluster::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  Cpptraj::Cluster::Metric_Data::DsArray cluster_dataset;
  DataSet_Coords* coords = 0;
  // First check for data
  std::string dataSetname = analyzeArgs.GetStringKey("data");
  if (!dataSetname.empty()) {
    ArgList dsnames(dataSetname, ",");
    DataSetList inputDsets;
    for (ArgList::const_iterator name = dsnames.begin(); name != dsnames.end(); ++name) {
      DataSetList tempDSL = setup.DSL().GetMultipleSets( *name );
      if (tempDSL.empty()) {
        mprinterr("Error: %s did not correspond to any data sets.\n", dataSetname.c_str());
        return Analysis::ERR;
      }
      inputDsets += tempDSL;
    }
    for (DataSetList::const_iterator ds = inputDsets.begin(); ds != inputDsets.end(); ++ds) {
      // Clustering only allowed on 1D data sets.
      if ( (*ds)->Group() != DataSet::SCALAR_1D ) {
        mprinterr("Error: Clustering only allowed on 1D scalar data sets, %s is %zuD.\n",
                  (*ds)->legend(), (*ds)->Ndim());
        return Analysis::ERR;
      }
      cluster_dataset.push_back( *ds );
    }

    if (control_.SetupForDataSets(cluster_dataset, analyzeArgs, setup.DSL(), setup.DFL(), debugIn))
      return Analysis::ERR;
  } else {
    // Attempt to get coords dataset from datasetlist
    std::string setname = analyzeArgs.GetStringKey("crdset");
    coords = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
    if (coords == 0) {
      mprinterr("Error: Could not locate COORDS set corresponding to %s\n",
                setname.c_str());
      return Analysis::ERR;
    }

    if (control_.SetupForCoordsDataSet(coords, analyzeArgs, setup.DSL(), setup.DFL(), debugIn))
      return Analysis::ERR;
  }

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

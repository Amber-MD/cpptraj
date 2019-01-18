#include "Control.h"
#include "../CpptrajStdio.h"
#include "../DataSet_Coords.h"

void Cpptraj::Cluster::Control::Help() {
  mprintf("[crdset <COORDS set>]\n");
}

int Cpptraj::Cluster::Control::SetupForDataSets(DataSetList const& DSL, ArgList& analyzeArgs,
                                                int verbose)
{
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringKey("crdset");
  if (!setname.empty()) {
    DataSet_Coords* coords = (DataSet_Coords*)DSL.FindCoordsSet( setname );
    if (coords == 0) {
      mprinterr("Error: Could not locate COORDS set corresponding to %s\n",
                setname.c_str());
      return Analysis::ERR;
    }
  }


  return 0;
}

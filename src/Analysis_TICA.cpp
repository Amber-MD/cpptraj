#include "Analysis_TICA.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
Analysis_TICA::Analysis_TICA() :
  TgtTraj_(0),
  lag_(0)
{}

// Analysis_TICA::Help()
void Analysis_TICA::Help() const {

}

// Analysis_TICA::Setup()
Analysis::RetType Analysis_TICA::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringKey("crdset");
  TgtTraj_ = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
  if (TgtTraj_ == 0) {
    mprinterr("Error: Could not locate COORDS set corresponding to %s\n",
              setname.c_str());
    Help();
    return Analysis::ERR;
  }
  // Other keywords
  lag_ = analyzeArgs.getKeyInt("lag", 1);

  mprintf("    TICA: Time independent correlation analysis.\n");
  mprintf("\tUsing coordinates from set '%s'\n", TgtTraj_->legend());
  mprintf("\tTime lag: %i frames.\n", lag_);

  return Analysis::OK;
}

// Analysis_TICA::Analyze()
Analysis::RetType Analysis_TICA::Analyze() {

}

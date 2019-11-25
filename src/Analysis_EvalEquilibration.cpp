#include "Analysis_EvalEquilibration.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"
#include "RPNcalc.h"

Analysis_EvalEquilibration::Analysis_EvalEquilibration() :
  Analysis(HIDDEN),
  setIn_(0),
  debug_(0)
{}

// Analysis_EvalEquilibration::Help()
void Analysis_EvalEquilibration::Help() const {
  mprintf("\n");
}

// Analysis_EvalEquilibration::Setup()
Analysis::RetType Analysis_EvalEquilibration::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  debug_ = debugIn;
  std::string setInName = analyzeArgs.GetStringKey("set");
  if (setInName.empty()) {
    mprinterr("Error: Must specify input data set.\n");
    return Analysis::ERR;
  }
  setIn_ = setup.DSL().GetDataSet( setInName );
  if (setIn_ == 0) {
    mprinterr("Error: '%s' matches no data sets.\n", setInName.c_str());
    return Analysis::ERR;
  }
  if (setIn_->Group() != DataSet::SCALAR_1D) {
    mprinterr("Error: '%s' is not a 1D scalar set.\n", setIn_->legend());
    return Analysis::ERR;
  }
  dsname_ = analyzeArgs.GetStringKey("name");
  if (dsname_.empty())
    dsname_ = setup.DSL().GenerateDefaultName("EvalEquil");

  mprintf("    EVALEQUILIBRATION: Evaluate equilibration of set '%s'\n", setIn_->legend());
  mprintf("\tOutput set name: %s\n", dsname_.c_str());

  return Analysis::OK;
}

// Analysis_EvalEquilibration::Analyze()
Analysis::RetType Analysis_EvalEquilibration::Analyze() {
  // First do a linear fit.
  if (setIn_->Size() < 2) {
    mprinterr("Error: Not enough data in '%s' to evaluate.\n", setIn_->legend());
    return Analysis::ERR;
  }
  double slope, intercept, correl;
  DataSet_1D const& DS = static_cast<DataSet_1D const&>( *setIn_ );
  CpptrajFile statsout;
  statsout.OpenWrite("");
  int err = DS.LinearRegression( slope, intercept, correl, &statsout );

  if (err != 0) {
    mprinterr("Error: Could not perform linear regression fit.\n");
    return Analysis::ERR;
  }

  // Determine relaxation direction
  int relaxationDir = 0;
  if (slope < 0)
    relaxationDir = -1;
  else if (slope > 0)
    relaxationDir = 1;

  // Special case: if slope was exactly 0 (should be rare). Consider this
  // equilibrated.
  

  RPNcalc calc;
  calc.SetDebug(debug_);

  //if (calc.ProcessExpression( dsname_ + "
  return Analysis::OK;
}

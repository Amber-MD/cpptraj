#include "Analysis_EvalEquilibration.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"
#include "RPNcalc.h"

Analysis_EvalEquilibration::Analysis_EvalEquilibration() :
  Analysis(HIDDEN),
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

  dsname_ = analyzeArgs.GetStringKey("name");
  if (dsname_.empty())
    dsname_ = setup.DSL().GenerateDefaultName("EvalEquil");

  if (inputSets_.AddSetsFromArgs( analyzeArgs, setup.DSL() ))
    return Analysis::ERR;

  mprintf("    EVALEQUILIBRATION: Evaluate equilibration of %zu sets.\n", inputSets_.size());
  mprintf("\tOutput set name: %s\n", dsname_.c_str());

  return Analysis::OK;
}

// Analysis_EvalEquilibration::Analyze()
Analysis::RetType Analysis_EvalEquilibration::Analyze() {
  CpptrajFile statsout;
  statsout.OpenWrite("");

  RPNcalc calc;
  calc.SetDebug(debug_);

  for (Array1D::const_iterator it = inputSets_.begin(); it != inputSets_.end(); ++it)
  {
    DataSet_1D const& DS = static_cast<DataSet_1D const&>( *(*it) );
    // First do a linear fit.
    if (DS.Size() < 2) {
      mprintf("Warning: Not enough data in '%s' to evaluate.\n", DS.legend());
      continue;
    }
    double slope, intercept, correl;
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

    //if (calc.ProcessExpression( dsname_ + "
  }
  return Analysis::OK;
}

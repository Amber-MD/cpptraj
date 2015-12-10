#include "Analysis_Multicurve.h"
#include "Analysis_CurveFit.h"
#include "CpptrajStdio.h"

void Analysis_Multicurve::Help() const {
  mprintf("\tset <dset> [set <dset> ...]\n");
  Analysis_CurveFit::HelpText();
}


Analysis::RetType Analysis_Multicurve::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  masterDSL_ = setup.DSL_Ptr();
  masterDFL_ = setup.DFL_Ptr();
  debug_ = debugIn;
  // Parse all 'set' arguments.
  std::string set_arg = analyzeArgs.GetStringKey("set");
  while (!set_arg.empty()) {
    inputDsets_.AddDataSets( setup.DSL().GetMultipleSets( set_arg ) );
    set_arg = analyzeArgs.GetStringKey("set");
  }
  if (inputDsets_.empty()) {
    mprinterr("Error: No data sets specified with 'set'\n");
    return Analysis::ERR;
  }
  args_ = analyzeArgs.RemainingArgs();

  mprintf("    MULTICURVE: Performing curve fitting on %zu sets.\n", inputDsets_.size());
  mprintf("\tUsing args: [%s]\n", args_.ArgLine());
  return Analysis::OK;
}

Analysis::RetType Analysis_Multicurve::Analyze() {
  int err = 0;
  for (Array1D::const_iterator set = inputDsets_.begin(); set != inputDsets_.end(); ++set) {
    ArgList argIn = args_;
    Analysis_CurveFit fit( (DataSet*)*set, set - inputDsets_.begin(), argIn,
                           *masterDSL_, *masterDFL_, debug_ );
    if (fit.Analyze()) ++err;
    mprintf("\n");
  }
  if (err > 0) return Analysis::ERR;
  return Analysis::OK;
}

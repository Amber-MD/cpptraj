#include "Analysis_ConstantPHStats.h"
#include "CpptrajStdio.h"

// Analysis_ConstantPHStats::Help()
void Analysis_ConstantPHStats::Help() const {

}

// Analysis_ConstantPHStats::Setup()
Analysis::RetType Analysis_ConstantPHStats::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  debug_ = debugIn;
  // Get DataSets
  std::string dsarg = analyzeArgs.GetStringNext();
  while (!dsarg.empty()) {
    inputSets_ += setup.DSL().GetMultipleSets( dsarg );
    dsarg = analyzeArgs.GetStringNext();
  }

  mprintf("    CONSTANT PH STATS:\n");
  inputSets_.List();
  return Analysis::OK;
}

// Analysis_ConstantPHStats::Analyze()
Analysis::RetType Analysis_ConstantPHStats::Analyze() {
  return Analysis::ERR;
}

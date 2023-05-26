#include "Analysis_CalcDiffusion.h"
#include "CpptrajStdio.h"

// Analysis_CalcDiffusion::Help()
void Analysis_CalcDiffusion::Help() const {
  mprintf("\t[crdset <coords set>] [maxlag <maxlag>] [<mask>]\n");
}

// Analysis_CalcDiffusion::Setup()
Analysis::RetType Analysis_CalcDiffusion::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringKey("crdset");
  TgtTraj_ = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
  if (TgtTraj_ == 0) {
    mprinterr("Error: Could not locate COORDS set corresponding to '%s'\n",
              setname.c_str());
    Help();
    return Analysis::ERR;
  }
  maxlag_ = analyzeArgs.getKeyInt("maxlag", -1);
  // Mask
  if (mask1_.SetMaskString( analyzeArgs.GetMaskNext() )) {
    mprinterr("Error: Could not set mask string.\n");
    return Analysis::ERR;
  }

  mprintf("    CALCDIFFUSION: Calculating diffusion from COORDS set '%s'\n", TgtTraj_->legend());
  if (maxlag_ > 0)
    mprintf("\tMaximum lag is %i frames.\n", maxlag_);
  mprintf("\tUsing atoms selected by mask '%s'\n", mask1_.MaskString());

  return Analysis::OK;
}

// Analysis_CalcDiffusion::Analyze()
Analysis::RetType Analysis_CalcDiffusion::Analyze() {

}

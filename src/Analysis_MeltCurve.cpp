#include "Analysis_MeltCurve.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"

void Analysis_MeltCurve::Help() const {
  mprintf("\t<dset0> [<dset1> ...] [out <outfile>] [name <outsetname>] cut <cut>\n"
          "  Calculate melting curve from input data sets assuming a simple 2-state\n"
          "  transition model, using data below <cut>as 'folded' and data above <cut>\n"
          "  as 'unfolded'\n");
}

// Analysis_MeltCurve::Setup()
Analysis::RetType Analysis_MeltCurve::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  std::string setname = analyzeArgs.GetStringKey("name");
  DataFile* outfile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  cut_ = analyzeArgs.getKeyDouble("cut", -1.0);
  if (cut_ < 0.0) {
    mprinterr("Error: meltcurve: 'cut <cut>' must be specified and > 0.0\n");
    return Analysis::ERR;
  }
  // Select datasets from remaining args
  if (input_dsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), setup.DSL() )) {
    mprinterr("Error: meltcurve: could not add data sets.\n");
    return Analysis::ERR;
  }

  // Set up output dataset
  mcurve_ = setup.DSL().AddSet(DataSet::DOUBLE, setname, "Melt");
  if (mcurve_ == 0) return Analysis::ERR;
  // Add dataset to datafile
  if (outfile != 0) outfile->AddDataSet( mcurve_ );
  
  mprintf("    MELTCURVE: Calculating melting curve from %zu data sets.\n",
          input_dsets_.size());
  mprintf("\tCut= %f", cut_);
  if (!setname.empty())
    mprintf("  Output set name: %s", setname.c_str());
  if (outfile != 0)
    mprintf("\tOutfile name: %s", outfile->DataFilename().base());
  mprintf("\n");
  for (Array1D::const_iterator set = input_dsets_.begin(); set != input_dsets_.end(); ++set)
    mprintf("\t%s\n", (*set)->legend());
  return Analysis::OK;
}

// Analysis_MeltCurve::Analyze()
Analysis::RetType Analysis_MeltCurve::Analyze() {
  int idx = 0;
  for (Array1D::const_iterator DS = input_dsets_.begin();
                               DS != input_dsets_.end(); ++DS)
  {
    if ( (*DS)->Size() < 1)
      mprintf("Warning: Set [%i] \"%s\" has no data.\n", idx, (*DS)->legend());
    else {
      int n_folded = 0;
      for (unsigned int i = 0; i < (*DS)->Size(); i++) {
        if ((*DS)->Dval(i) < cut_)
          ++n_folded;
      }
      double frac = (double)n_folded / (double)(*DS)->Size();
      mcurve_->Add(idx, &frac);
    }
    ++idx;
  }
  return Analysis::OK;
}

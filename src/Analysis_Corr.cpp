#include "Analysis_Corr.h"
#include "CpptrajStdio.h"
#include "DS_Math.h"

// CONSTRUCTOR
Analysis_Corr::Analysis_Corr() :
  D1_(0),
  D2_(0),
  lagmax_(0),
  usefft_(true),
  calc_covar_(true)
{}

void Analysis_Corr::Help() {
  mprintf("\tout <outfilename> <Dataset1> [<Dataset2>] [name <name>]\n");
  mprintf("\t[lagmax <lag>] [nocovar] [direct]\n");
  mprintf("\tCalculate auto-correlation for <Dataset1>, or cross -correlation\n");
  mprintf("\tbetween <Dataset1> and <Dataset2>\n");
}

// Analysis_Corr::Setup()
Analysis::RetType Analysis_Corr::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  const char* calctype;
  // Keywords
  lagmax_ = analyzeArgs.getKeyInt("lagmax",-1);
  usefft_ = !analyzeArgs.hasKey("direct");
  calc_covar_ = !analyzeArgs.hasKey("nocovar");
  DataFile* outfile = DFLin->AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  if (outfile == 0) {
    mprinterr("Error: Corr: No output filename specified ('out' <filename>).\n");
    return Analysis::ERR;
  }
 
  // DataSet names
  std::string D1name = analyzeArgs.GetStringNext();
  if (D1name.empty()) {
    mprinterr("Error: Corr: Must specify at least 1 dataset name.\n");
    return Analysis::ERR;
  }
  std::string D2name = analyzeArgs.GetStringNext();
  // Get DataSet(s)
  D1_ = datasetlist->GetDataSet(D1name);
  if (D1_==0) {
    mprinterr("Error: Corr: Could not get dataset named %s\n",D1name.c_str());
    return Analysis::ERR;
  }
  if (!D2name.empty())
    D2_ = datasetlist->GetDataSet(D2name);
  else {
    D2_ = D1_;
    D2name = D1name;
  }
  if (D2_==0) {
    mprinterr("Error: Corr: Could not get dataset named %s\n",D2name.c_str());
    return Analysis::ERR;
  }

  // TODO: Check DataSet type
  std::string dataset_name = analyzeArgs.GetStringKey("name");
  if (dataset_name.empty())
    dataset_name = datasetlist->GenerateDefaultName( "Corr" );

  // Setup output dataset
  std::string corrname = D1_->Legend() + "-" + D2_->Legend();
  Ct_ = datasetlist->AddSetAspect( DataSet::DOUBLE, dataset_name, corrname );
  if (Ct_ == 0) return Analysis::ERR;
  outfile->AddSet( Ct_ );

  if (calc_covar_)
    calctype = "covariance";
  else
    calctype = "correlation";

  if (D1name == D2name)
    mprintf("    CORR: auto-%s of set %s", calctype, D1name.c_str());
  else
    mprintf("    CORR: %s between set %s and set %s", calctype, D1name.c_str(), D2name.c_str());
  if (lagmax_!=-1) 
    mprintf(", max lag %i",lagmax_);
  mprintf("\n\tOutput to %s\n",outfile->DataFilename().base());
  if (usefft_)
    mprintf("\tUsing FFT to calculate %s.\n", calctype);
  else
    mprintf("\tUsing direct method to calculate %s.\n", calctype);

  return Analysis::OK;
}

// Analysis_Corr::Analyze()
Analysis::RetType Analysis_Corr::Analyze() {
  // Check that D1 and D2 have same # data points.
  int Nelements = D1_->Size(); 
  if (Nelements != D2_->Size()) {
    mprinterr("Error: Corr: # elements in dataset %s (%i) not equal to\n",
              D1_->Legend().c_str(), Nelements);
    mprinterr("             # elements in dataset %s (%i)\n",
              D2_->Legend().c_str(), D2_->Size());
    return Analysis::ERR;
  }
  if (lagmax_==-1) lagmax_ = Nelements;

  mprintf("    CORR: %i elements, max lag %i\n",Nelements,lagmax_);

  DS_Math::CrossCorr(*D1_, *D2_, *Ct_, lagmax_, calc_covar_, usefft_ );

  mprintf("    CORRELATION COEFFICIENT %6s to %6s IS %10.4f\n",
          D1_->Legend().c_str(), D2_->Legend().c_str(), DS_Math::CorrCoeff( *D1_, *D2_ ) );

  return Analysis::OK;
}

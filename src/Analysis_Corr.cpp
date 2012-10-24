#include "Analysis_Corr.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Analysis_Corr::Analysis_Corr() :
  D1_(NULL),
  D2_(NULL),
  lagmax_(0),
  usefft_(true),
  calc_covar_(true)
{}

void Analysis_Corr::Help() {
  mprintf("corr out <outfilename> <Dataset1> [<Dataset2>] [lagmax <lagmax>]\n");
}

// Analysis_Corr::Setup()
Analysis::RetType Analysis_Corr::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, int debugIn)
{
  const char* calctype;
  // Keywords
  lagmax_ = analyzeArgs.getKeyInt("lagmax",-1);
  usefft_ = !analyzeArgs.hasKey("direct");
  outfilename_ = analyzeArgs.GetStringKey("out");
  if (outfilename_.empty()) {
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
  if (D2_==NULL) {
    mprinterr("Error: Corr: Could not get dataset named %s\n",D2name.c_str());
    return Analysis::ERR;
  }

  // TODO: Check DataSet type
  std::string dataset_name = analyzeArgs.GetStringNext();
  if (dataset_name.empty())
    dataset_name = datasetlist->GenerateDefaultName( "Corr" );

  // Setup output dataset
  std::string corrname = D1_->Legend() + "-" + D2_->Legend();
  Ct_ = datasetlist->AddSetAspect( DataSet::DOUBLE, dataset_name, corrname );
  if (Ct_ == NULL) return Analysis::ERR;

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
  mprintf("\n\tOutput to %s\n",outfilename_.c_str());
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

  D1_->CrossCorr( *D2_, *Ct_, lagmax_, calc_covar_, usefft_ );

  mprintf("    CORRELATION COEFFICIENT %6s to %6s IS %10.4f\n",
          D1_->Legend().c_str(), D2_->Legend().c_str(), D1_->Corr( *D2_ ) );

  return Analysis::OK;
}

// Analysis_Corr::Print()
void Analysis_Corr::Print(DataFileList *dfl) {
  dfl->AddSetToFile(outfilename_, Ct_);
}


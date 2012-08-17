#include <cmath>
#include "Analysis_Corr.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Analysis_Corr::Analysis_Corr() :
  D1_(NULL),
  D2_(NULL),
  lagmax_(0)
{}

// Analysis_Corr::Setup()
/** Expected call: corr out <outfilename> <Dataset1> [<Dataset2>] [lagmax <lagmax>]
  */
int Analysis_Corr::Setup(DataSetList *datasetlist) {
  // If command was 'analyze correlationcoe' instead of 'corr' make sure
  // first two args are marked.
  if (analyzeArgs_[0] == "analyze") {
    analyzeArgs_.MarkArg(0);
    analyzeArgs_.MarkArg(1);
  }
  // Keywords
  lagmax_ = analyzeArgs_.getKeyInt("lagmax",-1);
  outfilename_ = analyzeArgs_.GetStringKey("out");
  if (outfilename_.empty()) {
    mprinterr("Error: Corr: No output filename specified ('out' <filename>).\n");
    return 1;
  }
 
  // DataSet names
  ArgList::ConstArg D1name = analyzeArgs_.getNextString();
  if (D1name==NULL) {
    mprinterr("Error: Corr: Must specify at least 1 dataset name.\n");
    return 1;
  }
  ArgList::ConstArg D2name = analyzeArgs_.getNextString();
  // Get DataSet(s)
  D1_ = datasetlist->Get(D1name);
  if (D1_==NULL) {
    mprinterr("Error: Corr: Could not get dataset named %s\n",D1name);
    return 1;
  }
  if (D2name!=NULL)
    D2_ = datasetlist->Get(D2name);
  else {
    D2_ = D1_;
    D2name = D1name;
  }
  if (D2_==NULL) {
    mprinterr("Error: Corr: Could not get dataset named %s\n",D2name);
    return 1;
  }

  // TODO: Check DataSet type

  // Setup output dataset
  std::string corrname = D1_->Legend() + "-" + D2_->Legend();
  Ct_ = datasetlist->AddSetAspect( DataSet::DOUBLE, "Corr", corrname );
  if (Ct_ == NULL) return 1;

  if (D1name == D2name)
    mprintf("    CORR: Auto-correlation of set %s", D1name);
  else
    mprintf("    CORR: Correlation between set %s and set %s",D1name,D2name);
  if (lagmax_!=-1) 
    mprintf(", max lag %i",lagmax_);
  mprintf("\n\tOutput to %s\n",outfilename_.c_str());

  return 0;
}

// Analysis_Corr::Analyze()
int Analysis_Corr::Analyze() {
  // Check that D1 and D2 have same # data points.
  int Nelements = D1_->Size(); 
  if (Nelements != D2_->Size()) {
    mprinterr("Error: Corr: # elements in dataset %s (%i) not equal to\n",D1_->c_str(),Nelements);
    mprinterr("             # elements in dataset %s (%i)\n",D2_->c_str(), D2_->Size());
    return 1;
  }
  if (lagmax_==-1) lagmax_ = Nelements;

  mprintf("    CORR: %i elements, max lag %i\n",Nelements,lagmax_);

  double corr_coeff = D1_->CrossCorr( *D2_, *Ct_, lagmax_, true );

  mprintf("    CORRELATION COEFFICIENT %6s to %6s IS %10.4f\n",
          D1_->c_str(), D2_->c_str(), corr_coeff );

  return 0;
}

// Analysis_Corr::Print()
void Analysis_Corr::Print(DataFileList *dfl) {
  dfl->AddSetToFile(outfilename_, Ct_);
}


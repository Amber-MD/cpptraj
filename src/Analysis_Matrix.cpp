#include "Analysis_Matrix.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Analysis_Matrix::Analysis_Matrix() :
  minfo_(0)
{}

// Command: analyze matrix
int Analysis_Matrix::Setup(DataSetList* DSLin) {
  // Ensure first 2 args (should be 'analyze' 'matrix') are marked.
  analyzeArgs_.MarkArg(0);
  analyzeArgs_.MarkArg(1);
  // Get matrix name
  std::string mname = analyzeArgs_.GetStringNext();
  if (mname.empty()) {
    mprinterr("Error: analyze matrix: missing matrix name (first argument).\n");
    return 1;
  }
  // Find matrix in DataSetList.
  DSLin->VectorBegin();
  while ( (minfo_ = (MatrixType*)DSLin->NextMatrix()) != 0 ) {
    if (mname == minfo_->Name() ) break;
  }
  if (minfo_ == 0) {
    mprinterr("Error: analyze matrix: Could not find matrix named %s\n",mname.c_str());
    return 1;
  }

  outfilename_ = analyzeArgs_.GetStringKey("out");
  orderparamfile_ = analyzeArgs_.GetStringKey("orderparamfile");
  if (analyzeArgs_.hasKey("thermo"))
    thermo_ = THERMO;
  else if (analyzeArgs_.hasKey("order"))
    thermo_ = ORDER;
  else
    thermo_ = OFF;
  if (thermo_==THERMO && minfo_->Type()!=MatrixType::MATRIX_MWCOVAR) {
    mprinterr("Error: analyze matrix: parameter 'thermo' only works for\n");
    mprinterr("       mass-weighted covariance matrix ('mwcovar').\n");
    return 1;
  }
  if (thermo_==ORDER && minfo_->Type()!=MatrixType::MATRIX_IRED) {
    mprinterr("Error: analyze matrix: parameter 'order' only works for\n");
    mprinterr("       IRED matrices.\n");
    return 1;
  }
  
  return 0;
}


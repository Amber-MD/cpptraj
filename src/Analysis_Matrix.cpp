#include "Analysis_Matrix.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Analysis_Matrix::Analysis_Matrix() :
  minfo_(0),
  thermo_(OFF),
  nevec_(0),
  reduce_(false),
  modinfo_(0)
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
  // Filenames
  outfilename_ = analyzeArgs_.GetStringKey("out");
  orderparamfile_ = analyzeArgs_.GetStringKey("orderparamfile");
  // Thermo/order flag
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
  // Number of eigenvectors; allow "0" only in case of 'thermo'
  nevec_ = analyzeArgs_.getKeyInt("vecs",0);
  if (thermo_ == THERMO) {
    if (nevec_ < 0) nevec_ = 0;
  } else if (nevec_ <= 0) {
    nevec_ = 1;
  }
  // Reduce flag
  reduce_ = analyzeArgs_.hasKey("reduce");
  // Set up ModesInfo
  std::string modesname = analyzeArgs_.GetStringKey("name");
  // TODO: Pass in MatrixType Nelt
  modinfo_ = new ModesInfo( (ModesInfo::modesType)minfo_->Type(),
                            minfo_->Mask1Tot(), ModesInfo::MS_STACK,
                            modesname );
  // Unlike PTRAJ, always put ModesInfo into the DataSetList; this way
  // DataSetList always handles the destruction.
  // TODO: Check for name conflicts.
  DSLin->AddDataSet( (DataSet*)modinfo_ );

  // Print Status
  mprintf("    ANALYZE MATRIX: Analyzing matrix %s",minfo_->c_str());
  if (!outfilename_.empty())
    mprintf(" and dumping results to %s\n", outfilename_.c_str());
  else
    mprintf(" and printing results to STDOUT\n");
  mprintf("      Calculating %i eigenvectors", nevec_);
  if (thermo_==THERMO)
    mprintf(" and thermodynamic data");
  else if (thermo_==ORDER)
    mprintf(" and order parameters");
  mprintf("\n");
  if (nevec_>0 && reduce_)
    mprintf("      Eigenvectors will be reduced\n");
  if (!modesname.empty())
    mprintf("      Storing modes with name: %s\n",modesname.c_str());
  if (!orderparamfile_.empty())
    mprintf("      Order parameters will be written to %s\n",orderparamfile_.c_str());
  return 0;
}


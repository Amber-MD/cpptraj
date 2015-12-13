#include "Analysis_CrossCorr.h"
#include "DataSet_MatrixFlt.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString
#include "DataSet_1D.h"

// CONSTRUCTOR
Analysis_CrossCorr::Analysis_CrossCorr() : outfile_(0), matrix_(0) {}

void Analysis_CrossCorr::Help() const {
  mprintf("\t[name <dsetname>] <dsetarg0> [<dsetarg1> ...] [out <filename>]\n"
          "  Calculate matrix of Pearson product-moment correlation\n"
          "  coefficients between selected data sets.\n");
}

// Analysis_CrossCorr::Setup()
Analysis::RetType Analysis_CrossCorr::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  std::string setname = analyzeArgs.GetStringKey("name");
  outfile_ = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  // Select datasets from remaining args
  if (input_dsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), setup.DSL() )) {
    mprinterr("Error: Could not add data sets.\n");
    return Analysis::ERR;
  }
  if (input_dsets_.size() < 2) {
    mprinterr("Error: At least 2 data sets are required.\n");
    return Analysis::ERR;
  }
  // Setup output dataset
  matrix_ = setup.DSL().AddSet( DataSet::MATRIX_FLT, setname, "crosscorr" );
  if (outfile_ != 0) {
    matrix_->SetDim(Dimension::X, Dimension(1.0, 1.0, "DataSets"));
    outfile_->AddDataSet( matrix_ );
  }
  
  mprintf("    CROSSCORR: Calculating correlation between %zu data sets:\n", input_dsets_.size());
  for (Array1D::const_iterator ds = input_dsets_.begin(); ds != input_dsets_.end(); ++ds)
    mprintf("\t'%s'\n", (*ds)->legend());
  mprintf("\tOutput set name: %s\n", matrix_->Meta().Name().c_str() );
  if ( outfile_ != 0 )
    mprintf("\tOutfile name: %s\n", outfile_->DataFilename().full());

  return Analysis::OK;
}

// Analysis_CrossCorr::Analyze()
Analysis::RetType Analysis_CrossCorr::Analyze() {
  DataSet_MatrixFlt& tmatrix = static_cast<DataSet_MatrixFlt&>( *matrix_ );
  if (tmatrix.AllocateTriangle(input_dsets_.size())) return Analysis::ERR;

  mprintf("\tDataSet Legend:\n");
  std::string Ylabels("\"");
  for (Array1D::const_iterator ds = input_dsets_.begin(); ds != input_dsets_.end(); ++ds) {
    int idx = (int)(ds - input_dsets_.begin() + 1);
    mprintf("\t\t%8i: %s\n", idx, (*ds)->legend());
    //Xlabels_ += (dsets_[i]->Legend() + ",");
    Ylabels += (integerToString(idx) + ":" + (*ds)->Meta().Legend() + ",");
  }
  Ylabels += "\"";
  for (Array1D::const_iterator ds0 = input_dsets_.begin(); ds0 != input_dsets_.end(); ++ds0) {
    DataSet_1D const& set0 = static_cast<DataSet_1D const&>( *(*ds0) );
    for (Array1D::const_iterator ds1 = ds0 + 1; ds1 != input_dsets_.end(); ++ds1) {
      float corr = (float)set0.CorrCoeff( *(*ds1) );
      //mprinterr("DBG:\tCross corr between %i (%s) and %i (%s)\n",
      //          i, dsets_[i]->legend(), j, dsets_[j]->legend());
      tmatrix.AddElement( corr );
    }
  }
  if (outfile_ != 0)
    outfile_->ProcessArgs("ylabels " + Ylabels);

  return Analysis::OK;
}

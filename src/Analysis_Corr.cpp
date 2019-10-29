#include "Analysis_Corr.h"
#include "CpptrajStdio.h"
#include "DataSet_Vector.h"

// CONSTRUCTOR
Analysis_Corr::Analysis_Corr() :
  D1_(0),
  D2_(0),
  Ct_(0),
  Coeff_(0),
  lagmax_(0),
  usefft_(true),
  calc_covar_(true)
{}

void Analysis_Corr::Help() const {
  mprintf("\tout <outfilename> <Dataset1> [<Dataset2>] [name <name>]\n"
          "\t[lagmax <lag>] [nocovar] [direct]\n"
          "  Calculate auto-correlation for <Dataset1>, or cross -correlation\n"
          "  between <Dataset1> and <Dataset2>\n");
}

// Analysis_Corr::Setup()
Analysis::RetType Analysis_Corr::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  const char* calctype;
  // Keywords
  lagmax_ = analyzeArgs.getKeyInt("lagmax",-1);
  usefft_ = !analyzeArgs.hasKey("direct");
  calc_covar_ = !analyzeArgs.hasKey("nocovar");
  DataFile* outfile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  if (outfile == 0) {
    mprinterr("Error: Corr: No output filename specified ('out' <filename>).\n");
    return Analysis::ERR;
  }
  // TODO: Check DataSet type
  std::string dataset_name = analyzeArgs.GetStringKey("name");

  // DataSet names
  std::string D1name = analyzeArgs.GetStringNext();
  if (D1name.empty()) {
    mprinterr("Error: Corr: Must specify at least 1 dataset name.\n");
    return Analysis::ERR;
  }
  std::string D2name = analyzeArgs.GetStringNext();
  // Get DataSet(s)
  D1_ = setup.DSL().GetDataSet(D1name);
  if (D1_==0) {
    mprinterr("Error: Corr: Could not get dataset named %s\n",D1name.c_str());
    return Analysis::ERR;
  }
  if (!D2name.empty())
    D2_ = setup.DSL().GetDataSet(D2name);
  else {
    D2_ = D1_;
    D2name = D1name;
  }
  if (D2_==0) {
    mprinterr("Error: Corr: Could not get dataset named %s\n",D2name.c_str());
    return Analysis::ERR;
  }
  if (D1_->Type() == DataSet::VECTOR && D2_->Type() != DataSet::VECTOR)
  {
    mprinterr("Error: Vector cross correlation requires 2 vector data sets.\n");
    return Analysis::ERR;
  }

  // Setup output datasets
  Ct_ = setup.DSL().AddSet( DataSet::DOUBLE, dataset_name, "Corr" );
  if (Ct_ == 0) return Analysis::ERR;
  if (dataset_name.empty()) {
    std::string corrname = "C(" + D1_->Meta().Legend();
    if (D2_ != D1_) corrname += ("-" + D2_->Meta().Legend());
    corrname += ")";
    Ct_->SetLegend( corrname );
  }
  outfile->AddDataSet( Ct_ );
  Coeff_ = 0;
  if (D1_->Type() != DataSet::VECTOR) {
    Coeff_ = setup.DSL().AddSet( DataSet::DOUBLE, MetaData(Ct_->Meta().Name(), "coeff") );
    if (Coeff_ == 0) return Analysis::ERR;
    Coeff_->Allocate( DataSet::SizeArray(1, 1) );
  }

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
  mprintf("\tCorrelation function data set: %s\n", Ct_->Meta().PrintName().c_str());
  if (Coeff_ != 0)
    mprintf("\tCorrelation coefficient data set: %s\n", Coeff_->Meta().PrintName().c_str());

  return Analysis::OK;
}

// Analysis_Corr::Analyze()
Analysis::RetType Analysis_Corr::Analyze() {
  // Check that D1 and D2 have same # data points.
  size_t Nelements = D1_->Size(); 
  if (Nelements != D2_->Size()) {
    mprinterr("Error: Corr: # elements in dataset %s (%zu) not equal to\n",
              D1_->legend(), Nelements);
    mprinterr("             # elements in dataset %s (%zu)\n",
              D2_->legend(), D2_->Size());
    return Analysis::ERR;
  }
  if (lagmax_==-1) lagmax_ = (int)Nelements;

  mprintf("    CORR: %zu elements, max lag %i\n",Nelements,lagmax_);

  Analysis::RetType err = Analysis::OK;
  if (D1_->Type() == DataSet::VECTOR) {
    DataSet_Vector const& set1 = static_cast<DataSet_Vector const&>( *D1_ );
    DataSet_Vector const& set2 = static_cast<DataSet_Vector const&>( *D2_ );
    set1.CalcVectorCorr(set2, *((DataSet_1D*)Ct_), lagmax_);
  } else {
    DataSet_1D const& set1 = static_cast<DataSet_1D const&>( *D1_ );
    DataSet_1D const& set2 = static_cast<DataSet_1D const&>( *D2_ );
    if (set1.CrossCorr( set2, *((DataSet_1D*)Ct_), lagmax_, calc_covar_, usefft_ ))
      err = Analysis::ERR;
    double coeff = set1.CorrCoeff( set2 );
    mprintf("    CORRELATION COEFFICIENT %s to %s IS %.4f\n",
            D1_->legend(), D2_->legend(), coeff );
    Coeff_->Add(0, &coeff);
  }

  return err;
}

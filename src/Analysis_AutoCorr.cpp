#include "Analysis_AutoCorr.h"
#include "CpptrajStdio.h"
#include "DataSet_Vector.h"

// CONSTRUCTOR
Analysis_AutoCorr::Analysis_AutoCorr() :
  lagmax_(-1),
  usefft_(true),
  calc_covar_(true)
{}

void Analysis_AutoCorr::Help() const {
  mprintf("\t[name <dsetname>] <dsetarg0> [<dsetarg1> ...] [out <filename>]\n"
          "\t[lagmax <lag>] [nocovar] [direct]\n"
          "  Calculate autocorrelation functions for selected data set(s)\n");
}

Analysis::RetType Analysis_AutoCorr::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  const char* calctype;

  std::string setname = analyzeArgs.GetStringKey("name");
  DataFile* outfile = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
  lagmax_ = analyzeArgs.getKeyInt("lagmax",-1);
  calc_covar_ = !analyzeArgs.hasKey("nocovar");
  usefft_ = !analyzeArgs.hasKey("direct");
  // Select datasets from remaining args
  dsets_.clear();
  ArgList dsetArgs = analyzeArgs.RemainingArgs();
  for (ArgList::const_iterator dsa = dsetArgs.begin(); dsa != dsetArgs.end(); ++dsa) {
    DataSetList setsIn = setup.DSL().GetMultipleSets( *dsa );
    for (DataSetList::const_iterator ds = setsIn.begin(); ds != setsIn.end(); ++ds) {
      if ( (*ds)->Group() != DataSet::SCALAR_1D && (*ds)->Type() != DataSet::VECTOR )
        mprintf("Warning: Set '%s' type not supported in AUTOCORR - skipping.\n",
                (*ds)->legend());
      else
        dsets_.push_back( *ds );
    }
  }
  if (dsets_.empty()) {
    mprinterr("Error: No data sets selected.\n");
    return Analysis::ERR;
  }
  // If setname is empty generate a default name
  if (setname.empty())
    setname = setup.DSL().GenerateDefaultName( "autocorr" );
  // Setup output datasets
  MetaData md( setname );
  for (unsigned int idx = 0; idx != dsets_.size(); idx++) {
    md.SetIdx( idx );
    DataSet* dsout = setup.DSL().AddSet( DataSet::DOUBLE, md );
    if (dsout==0) return Analysis::ERR;
    dsout->SetLegend( dsets_[idx]->Meta().Legend() );
    outputData_.push_back( dsout );
    // Add set to output file
    if (outfile != 0) outfile->AddDataSet( outputData_.back() );
  }
 
  if (calc_covar_)
    calctype = "covariance";
  else
    calctype = "correlation";
 
  mprintf("    AUTOCORR: Calculating auto-%s for %i data sets:\n\t", calctype, dsets_.size());
  for (unsigned int idx = 0; idx != dsets_.size(); ++idx)
    mprintf(" %s", dsets_[idx]->legend());
  mprintf("\n");
  if (lagmax_!=-1)
    mprintf("\tLag max= %i\n", lagmax_);
  if ( !setname.empty() )
    mprintf("\tSet name: %s\n", setname.c_str() );
  if ( outfile != 0 )
    mprintf("\tOutfile name: %s\n", outfile->DataFilename().base());
  if (usefft_)
    mprintf("\tUsing FFT to calculate %s.\n", calctype);
  else
    mprintf("\tUsing direct method to calculate %s.\n", calctype);

  return Analysis::OK;
}

Analysis::RetType Analysis_AutoCorr::Analyze() {
  for (unsigned int ids = 0; ids != dsets_.size(); ids++) {
    mprintf("\t\tCalculating AutoCorrelation for set %s\n", dsets_[ids]->legend());
    DataSet_1D& Ct = static_cast<DataSet_1D&>( *outputData_[ids] );
    if (dsets_[ids]->Type() == DataSet::VECTOR) {
      DataSet_Vector const& set = static_cast<DataSet_Vector const&>( *dsets_[ids] );
      set.CalcVectorCorr( set, Ct, lagmax_ );
    } else {
      DataSet_1D const& set = static_cast<DataSet_1D const&>( *dsets_[ids] );
      set.CrossCorr( set, Ct, lagmax_, calc_covar_, usefft_ );
    }
  }

  return Analysis::OK;
}

#include "Analysis_AutoCorr.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Analysis_AutoCorr::Analysis_AutoCorr() :
  lagmax_(-1),
  usefft_(true),
  calc_covar_(true)
{}

/** Usage: autocorr [name <dsetname>] <dsetarg0> [<dsetarg1> ...] out <filename>
  */
int Analysis_AutoCorr::Setup( DataSetList* datasetlist ) {
  const char* calctype;

  std::string setname_ = analyzeArgs_.GetStringKey("name");
  outfilename_ = analyzeArgs_.GetStringKey("out");
  lagmax_ = analyzeArgs_.getKeyInt("lagmax",-1);
  // Select datasets
  while (analyzeArgs_.ArgsRemain())
    dsets_ += datasetlist->GetMultipleSets( analyzeArgs_.GetStringNext() );
  if (dsets_.empty()) {
    mprinterr("Error: autocorr: No data sets selected.\n");
    return 1;
  }
  // If setname is empty generate a default name
  if (setname_.empty())
    setname_ = datasetlist->GenerateDefaultName( "autocorr" );
  // Setup output datasets
  int idx = 0;
  for (DataSetList::const_iterator DS = dsets_.begin(); DS != dsets_.end(); ++DS) {
    DataSet* dsout = datasetlist->AddSetIdx( DataSet::DOUBLE, setname_, idx++ );
    if (dsout==NULL) return 1;
    dsout->SetLegend( (*DS)->Legend() );
    outputData_.push_back( dsout );
  }
 
  if (calc_covar_)
    calctype = "covariance";
  else
    calctype = "correlation";
 
  mprintf("    AUTOCORR: Calculating auto-%s for %i data sets:\n", calctype, dsets_.size());
  dsets_.Info();
  if (lagmax_!=-1)
    mprintf("\tLag max= %i\n", lagmax_);
  if ( !setname_.empty() )
    mprintf("\tSet name: %s\n", setname_.c_str() );
  if ( !outfilename_.empty() )
    mprintf("\tOutfile name: %s\n", outfilename_.c_str());
  if (usefft_)
    mprintf("\tUsing FFT to calculate %s.\n", calctype);
  else
    mprintf("\tUsing direct method to calculate %s.\n", calctype);

  return 0;
}

int Analysis_AutoCorr::Analyze() {
  std::vector<DataSet*>::iterator dsout = outputData_.begin();
  for (DataSetList::const_iterator DS = dsets_.begin(); DS != dsets_.end(); ++DS)
  {
    mprintf("\t\tCalculating AutoCorrelation for set %s\n", (*DS)->Legend().c_str());
    (*DS)->CrossCorr( *(*DS), *(*dsout), lagmax_, calc_covar_, usefft_ );
    ++dsout;
  }

  return 0;
}

void Analysis_AutoCorr::Print( DataFileList* datafilelist ) {
  if (!outfilename_.empty()) {
    for (std::vector<DataSet*>::iterator dsout = outputData_.begin();
                                         dsout != outputData_.end(); ++dsout)
      datafilelist->Add( outfilename_.c_str(), *dsout );
    //DataFile* DF = datafilelist->GetDataFile( outfilename_.c_str());
    //if (DF != NULL) 
    //  DF->ProcessArgs("xlabel DataSets");
  }
}


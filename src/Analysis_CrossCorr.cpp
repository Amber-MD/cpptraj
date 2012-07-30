#include "Analysis_CrossCorr.h"
#include "TriangleMatrix.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Analysis_CrossCorr::Analysis_CrossCorr() {}

/** Usage: crosscorr [name <dsetname>] <dsetarg0> [<dsetarg1> ...] out <filename>
  */
int Analysis_CrossCorr::Setup( DataSetList* datasetlist ) {
  std::string setname_ = analyzeArgs_.GetStringKey("name");
  outfilename_ = analyzeArgs_.GetStringKey("out");
  // Select datasets
  dsets_ = datasetlist->GetMultipleSets( analyzeArgs_.GetStringNext() );
  if (dsets_.empty()) {
    mprinterr("Error: crosscorr: No data sets selected.\n");
    return 1;
  }
  // Setup output dataset
  matrix_ = datasetlist->AddSet( DataSet::TRIMATRIX, setname_, "crosscorr" );
  
  mprintf("    CROSSCORR: Calculating correlation between %i data sets:\n", dsets_.size());
  dsets_.Info();
  if ( !setname_.empty() )
    mprintf("\tSet name: %s\n", setname_.c_str() );
  if ( !outfilename_.empty() )
    mprintf("\tOutfile name: %s\n", outfilename_.c_str());

  return 0;
}

int Analysis_CrossCorr::Analyze() {
  TriangleMatrix* tmatrix = (TriangleMatrix*)matrix_;

  int Nsets = dsets_.size();
  mprintf("\tDataSet Legend:\n");
  for (int i = 0; i < Nsets; ++i)
    mprintf("\t\t%8i: %s\n", i+1, dsets_[i]->Legend().c_str());
  int Nsets1 = Nsets - 1;
  tmatrix->Setup(Nsets);
  for (int i = 0; i < Nsets1; ++i) {
    for (int j = i + 1; j < Nsets; ++j) {
      //mprinterr("DBG:\tCross corr between %i (%s) and %i (%s)\n",
      //          i, dsets_[i]->Legend().c_str(), j, dsets_[j]->Legend().c_str());
      double corr = dsets_[i]->Corr( *dsets_[j], NULL, 0 );
      tmatrix->AddElement( corr );
    }
  }

  return 0;
}

void Analysis_CrossCorr::Print( DataFileList* datafilelist ) {
  if (!outfilename_.empty())
    datafilelist->Add( outfilename_.c_str(), matrix_ );
}


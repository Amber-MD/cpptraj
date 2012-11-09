#include "Analysis_CrossCorr.h"
#include "TriangleMatrix.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString

// CONSTRUCTOR
Analysis_CrossCorr::Analysis_CrossCorr() {}

void Analysis_CrossCorr::Help() {
  mprintf("crosscorr [name <dsetname>] <dsetarg0> [<dsetarg1> ...] out <filename>\n");
}

Analysis::RetType Analysis_CrossCorr::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, int debugIn)
{
  std::string setname_ = analyzeArgs.GetStringKey("name");
  outfilename_ = analyzeArgs.GetStringKey("out");
  // Select datasets from remaining args
  ArgList dsetArgs = analyzeArgs.RemainingArgs();
  for (ArgList::const_iterator dsa = dsetArgs.begin(); dsa != dsetArgs.end(); ++dsa)
    dsets_ += datasetlist->GetMultipleSets( *dsa );
  if (dsets_.empty()) {
    mprinterr("Error: crosscorr: No data sets selected.\n");
    return Analysis::ERR;
  }
  // Setup output dataset
  matrix_ = datasetlist->AddSet( DataSet::TRIMATRIX, setname_, "crosscorr" );
  
  mprintf("    CROSSCORR: Calculating correlation between %i data sets:\n", dsets_.size());
  dsets_.List();
  if ( !setname_.empty() )
    mprintf("\tSet name: %s\n", setname_.c_str() );
  if ( !outfilename_.empty() )
    mprintf("\tOutfile name: %s\n", outfilename_.c_str());

  return Analysis::OK;
}

Analysis::RetType Analysis_CrossCorr::Analyze() {
  TriangleMatrix* tmatrix = (TriangleMatrix*)matrix_;

  int Nsets = dsets_.size();
  mprintf("\tDataSet Legend:\n");
  Ylabels_.assign("\"");
  for (int i = 0; i < Nsets; ++i) {
    mprintf("\t\t%8i: %s\n", i+1, dsets_[i]->Legend().c_str());
    //Xlabels_ += (dsets_[i]->Legend() + ",");
    Ylabels_ += (integerToString(i+1) + ":" + dsets_[i]->Legend() + ",");
  }
  Ylabels_ += "\"";
  int Nsets1 = Nsets - 1;
  tmatrix->Setup(Nsets);
  for (int i = 0; i < Nsets1; ++i) {
    for (int j = i + 1; j < Nsets; ++j) {
      //mprinterr("DBG:\tCross corr between %i (%s) and %i (%s)\n",
      //          i, dsets_[i]->Legend().c_str(), j, dsets_[j]->Legend().c_str());
      double corr = dsets_[i]->Corr( *dsets_[j] );
      tmatrix->AddElement( corr );
    }
  }

  return Analysis::OK;
}

void Analysis_CrossCorr::Print( DataFileList* datafilelist ) {
  if (!outfilename_.empty()) {
    DataFile* DF = datafilelist->AddSetToFile( outfilename_, matrix_ );
    if (DF != NULL) {
      //DF->ProcessArgs("xlabels " + Xlabels_);
      DF->ProcessArgs("xlabel DataSets ylabels " + Ylabels_);
    }
  }
}


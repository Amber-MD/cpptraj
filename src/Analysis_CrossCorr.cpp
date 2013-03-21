#include "Analysis_CrossCorr.h"
#include "TriangleMatrix.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString
#include "DS_Math.h"

// CONSTRUCTOR
Analysis_CrossCorr::Analysis_CrossCorr() : outfile_(0), matrix_(0) {}

void Analysis_CrossCorr::Help() {
  mprintf("\t[name <dsetname>] <dsetarg0> [<dsetarg1> ...] [out <filename>]\n");
  mprintf("\tCalculate matrix of Pearson product-moment correlation\n");
  mprintf("\tcoefficients between selected data sets.\n");
}

Analysis::RetType Analysis_CrossCorr::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  std::string setname = analyzeArgs.GetStringKey("name");
  outfile_ = DFLin->AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  // Select datasets from remaining args
  ArgList dsetArgs = analyzeArgs.RemainingArgs();
  for (ArgList::const_iterator dsa = dsetArgs.begin(); dsa != dsetArgs.end(); ++dsa)
    dsets_ += datasetlist->GetMultipleSets( *dsa );
  if (dsets_.empty()) {
    mprinterr("Error: crosscorr: No data sets selected.\n");
    return Analysis::ERR;
  }
  // Setup output dataset
  matrix_ = datasetlist->AddSet( DataSet::TRIMATRIX, setname, "crosscorr" );
  if (outfile_ != 0) {
    outfile_->AddSet( matrix_ );
    outfile_->ProcessArgs("xlabel DataSets");
  }
  
  mprintf("    CROSSCORR: Calculating correlation between %i data sets:\n", dsets_.size());
  dsets_.List();
  if ( !setname.empty() )
    mprintf("\tSet name: %s\n", setname.c_str() );
  if ( outfile_ != 0 )
    mprintf("\tOutfile name: %s\n", outfile_->DataFilename().base());

  return Analysis::OK;
}

Analysis::RetType Analysis_CrossCorr::Analyze() {
  TriangleMatrix* tmatrix = (TriangleMatrix*)matrix_;

  int Nsets = dsets_.size();
  mprintf("\tDataSet Legend:\n");
  std::string Ylabels("\"");
  for (int i = 0; i < Nsets; ++i) {
    mprintf("\t\t%8i: %s\n", i+1, dsets_[i]->Legend().c_str());
    //Xlabels_ += (dsets_[i]->Legend() + ",");
    Ylabels += (integerToString(i+1) + ":" + dsets_[i]->Legend() + ",");
  }
  Ylabels += "\"";
  int Nsets1 = Nsets - 1;
  tmatrix->Setup(Nsets);
  for (int i = 0; i < Nsets1; ++i) {
    for (int j = i + 1; j < Nsets; ++j) {
      //mprinterr("DBG:\tCross corr between %i (%s) and %i (%s)\n",
      //          i, dsets_[i]->Legend().c_str(), j, dsets_[j]->Legend().c_str());
      double corr = DS_Math::CorrCoeff( *dsets_[i], *dsets_[j] );
      tmatrix->AddElement( corr );
    }
  }
  if (outfile_ != 0)
    outfile_->ProcessArgs("ylabels " + Ylabels);

  return Analysis::OK;
}

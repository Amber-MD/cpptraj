#include "Analysis_Lifetime.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Analysis_Lifetime::Analysis_Lifetime() :
  windowSize_(1)
{}

/** Usage: lifetime out <filename> name <setname> <dsetarg0> [ <dsetarg1> ... ]
  */
int Analysis_Lifetime::Setup( DataSetList* datasetlist ) {
  // Get Keywords
  outfilename_ = analyzeArgs_.GetStringKey("out");
  std::string setname_ = analyzeArgs_.GetStringKey("name");
  if (setname_.empty()) {
    mprinterr("Error: lifetime: No setname given 'name <setname>'\n");
    return 1;
  }
  windowSize_ = analyzeArgs_.getKeyInt("window", 1);
  if (windowSize_ < 1) {
    mprinterr("Error: lifetime: window <windowsize>: Must be > 0 (%i)\n", windowSize_);
    return 1;
  }
  // Select datasets
  inputDsets_ = datasetlist->GetMultipleSets( analyzeArgs_.GetStringNext() );
  if (inputDsets_.empty()) {
    mprinterr("Error: lifetime: No data sets selected.\n");
    return 1;
  }

  // Create output datasets
  int didx = 0;
  for (DataSetList::const_iterator set = inputDsets_.begin(); set != inputDsets_.end(); ++set)
  {
    DataSet* outSet = datasetlist->AddSetIdx( DataSet::FLOAT, setname_, didx++ );
    if (outSet==NULL) {
      mprinterr("Error: lifetime: Could not allocate output set for %s\n", 
                (*set)->Legend().c_str());
      return 1;
    }
    outSet->SetLegend( (*set)->Legend() );
    outputDsets_.push_back( outSet );
  }

  mprintf("    LIFETIME: Calculating avg lifetime of data in %i sets:\n", 
          inputDsets_.size());
  inputDsets_.Info();
  if (!outfilename_.empty())
    mprintf("\tOutfile name: %s\n", outfilename_.c_str());
  mprintf("\tWindow size for averaging: %i\n", windowSize_);

  return 0;
}

int Analysis_Lifetime::Analyze() {
  std::vector<DataSet*>::iterator outSet = outputDsets_.begin();
  for (DataSetList::const_iterator inSet = inputDsets_.begin(); 
                                   inSet != inputDsets_.end(); ++inSet)
  {
    mprintf("\t\tCalculating lifetimes for set %s\n", (*inSet)->Legend().c_str());
    // Loop over all values in set.
    double sum = 0;
    int windowcount = 0;
    int frame = 0;
    int setSize = (*inSet)->Size();
    for (int i = 0; i < setSize; ++i) {
      sum += (*inSet)->Dval(i);
      ++windowcount;
      if (windowcount == windowSize_) {
        float fval = (float)sum / (float)windowcount;
        //(*outSet)->Add( frame++, &fval );
        //(*outSet)->Add( i, &fval );
        for (int j = frame; j < frame+windowcount; ++j)
          (*outSet)->Add( j, &fval );
        frame += windowcount;
        windowcount = 0;
        sum = 0;
      }
    }
    ++outSet;
  }
  return 0;
}

void Analysis_Lifetime::Print(DataFileList* datafilelist) {
  if (!outfilename_.empty()) {
    for (std::vector<DataSet*>::iterator outSet = outputDsets_.begin();
                                         outSet != outputDsets_.end(); ++outSet)
      datafilelist->Add( outfilename_.c_str(), *outSet );
  }
}


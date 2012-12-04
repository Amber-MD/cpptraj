#include "Analysis_RunningAvg.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Analysis_RunningAvg::Analysis_RunningAvg()
{}

void Analysis_RunningAvg::Help() {
  mprintf("runningavg <dset1> [<dset2> ...] [name <dsetname>] [out <filename>]\n");
}

Analysis::RetType Analysis_RunningAvg::Setup(ArgList& analyzeArgs, 
                                             DataSetList* datasetlist,
                                             TopologyList* PFLin,
                                             int debugIn)
{
  outfilename_ = analyzeArgs.GetStringKey("out");
  setname_ = analyzeArgs.GetStringKey("name");
  // The remaining arguments are the data sets to take running averages of
  ArgList dsetArgs = analyzeArgs.RemainingArgs();
  // Build the data set list
  for (ArgList::const_iterator dsa = dsetArgs.begin(); 
                dsa != dsetArgs.end(); ++dsa)
    dsets_ += datasetlist->GetMultipleSets( *dsa );
  
  if (dsets_.empty()) {
    mprinterr("Error: runningavg: No data sets selected.\n");
    return Analysis::ERR;
  }
  // If setname is empty, generate a default name
  if (setname_.empty())
    setname_ = datasetlist->GenerateDefaultName( "runningavg" );
  // Setup output datasets
  int idx = 0;
  for (DataSetList::const_iterator DS = dsets_.begin(); DS != dsets_.end(); ++DS) {
    DataSet* dsout = datasetlist->AddSetIdx( DataSet::DOUBLE, setname_, idx++);
    if (dsout == NULL)
      return Analysis::ERR;
    dsout->SetLegend( (*DS)->Legend() );
    outputData_.push_back( dsout );
  }

  mprintf("    RUNNINGAVG: Calculating the running average for %i data sets:\n",
          dsets_.size());
  dsets_.List();
  if ( !outfilename_.empty() )
    mprintf("\tOutfile name: %s\n", outfilename_.c_str());

  return Analysis::OK;
}

Analysis::RetType Analysis_RunningAvg::Analyze() {
  std::vector<DataSet*>::iterator dsout = outputData_.begin();
  for (DataSetList::const_iterator DS = dsets_.begin(); DS != dsets_.end(); DS++)
  {
    mprintf("\t\tCalculating Running Average for set %s\n", (*DS)->Legend().c_str());
    double running_sum = 0;
    for (int i = 0; i < (*DS)->Size(); i++) {
      running_sum += (*DS)->Dval(i);
      double avg = running_sum / (i+1);
      (*dsout)->Add(i, &avg);
    }
    dsout++;
  }

  return Analysis::OK;
}

void Analysis_RunningAvg::Print( DataFileList* datafilelist) {
  if (!outfilename_.empty()) {
    for (std::vector<DataSet*>::iterator dsout = outputData_.begin();
                                         dsout != outputData_.end(); dsout++)
      datafilelist->AddSetToFile( outfilename_, *dsout );
  }
}

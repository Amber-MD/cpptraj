#include "Analysis_RunningAvg.h"
#include "CpptrajStdio.h"

#define MIN(X, Y) ( ( (X) < (Y) ) ? (X) : (Y) )
#define MAX(X, Y) ( ( (X) < (Y) ) ? (Y) : (X) )

// CONSTRUCTOR
Analysis_RunningAvg::Analysis_RunningAvg() :
  cumulative_(false),
  window_(5)
{}

void Analysis_RunningAvg::Help() {
  mprintf("\t<dset1> [<dset2> ...] [name <dsetname>] [out <filename>]\n");
  mprintf("\t[ [cumulative] | [window <window>] ]\n");
  mprintf("\tCalculate running average of data in selected data set(s)\n");
}

Analysis::RetType Analysis_RunningAvg::Setup(ArgList& analyzeArgs, DataSetList* datasetlist, 
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  DataFile* outfile = DFLin->AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  std::string setname = analyzeArgs.GetStringKey("name");
  cumulative_ = analyzeArgs.hasKey("cumulative");
  window_ = analyzeArgs.getKeyDouble("window", 5);

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
  if (setname.empty())
    setname = datasetlist->GenerateDefaultName( "runningavg" );
  // Setup output datasets
  int idx = 0;
  for (DataSetList::const_iterator DS = dsets_.begin(); DS != dsets_.end(); ++DS) {
    DataSet* dsout = datasetlist->AddSetIdx( DataSet::DOUBLE, setname, idx++);
    if (dsout == 0)
      return Analysis::ERR;
    dsout->SetLegend( (*DS)->Legend() );
    outputData_.push_back( dsout );
    if (outfile != 0) outfile->AddSet( dsout );
  }

  if (cumulative_)
    mprintf("RUNNINGAVG: Calculating the cumulative running average for %i data sets",
            dsets_.size());
  else
    mprintf("RUNNINGAVG: Calculating the running average for %i data sets with a %d-element window",
            dsets_.size(), window_);
  dsets_.List();
  if ( outfile != 0 )
    mprintf("\tOutfile name: %s\n", outfile->DataFilename().base());

  return Analysis::OK;
}

Analysis::RetType Analysis_RunningAvg::Analyze() {
  std::vector<DataSet*>::iterator dsout = outputData_.begin();
  for (DataSetList::const_iterator DS = dsets_.begin(); DS != dsets_.end(); DS++)
  {
    if (cumulative_) {
      mprintf("\t\tCalculating Cumulative Running Average for set %s\n",
              (*DS)->Legend().c_str());
      double running_sum = 0;
      for (int i = 0; i < (*DS)->Size(); i++) {
        running_sum += (*DS)->Dval(i);
        double avg = running_sum / (i+1);
        (*dsout)->Add(i, &avg);
      }
    }else {
      mprintf("\t\tCalculating Running Average for set %s\n",
              (*DS)->Legend().c_str());
      for (int i = 0; i < (*DS)->Size(); i++) {
        int npts = 0;
        double sum = 0;
        for (int j = MAX(0, i-window_); j < MIN((*DS)->Size(), i+window_); j++) {
          sum += (*DS)->Dval(i);
          npts++;
        }
        double avg = sum / (double) npts;
        (*dsout)->Add(i, &avg);
      }
    }
    dsout++;
  }

  return Analysis::OK;
}

#include "Analysis_Average.h"
#include "CpptrajStdio.h"

void Analysis_Average::Help() {
  mprintf("\t<dset0> [<dset1> ...] [out <file>] [noheader]\n"
          "  Calculate the average, standard deviation, min, and max of given data sets.\n");
}

// Analysis_Average::Setup()
Analysis::RetType Analysis_Average::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  writeHeader_ = !analyzeArgs.hasKey("noheader");
  outfile_ = DFLin->AddCpptrajFile(analyzeArgs.GetStringKey("out"), "DataSet Average",
                                   DataFileList::TEXT, true);
  if (outfile_ == 0) return Analysis::ERR;
  // Select datasets from remaining args
  if (input_dsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), *datasetlist )) {
    mprinterr("Error: Could not add data sets.\n");
    return Analysis::ERR;
  }
  if (input_dsets_.empty()) {
    mprinterr("Error: No input data sets.\n");
    return Analysis::ERR;
  }

  mprintf("    AVERAGE: Calculating average of %i data sets.\n",
          input_dsets_.size());
  mprintf("\tWriting results to %s\n", outfile_->Filename().full());
  //for (Array1D::const_iterator set = input_dsets_.begin(); set != input_dsets_.end(); ++set)
  //  mprintf("\t%s\n", (*set)->legend());

  return Analysis::OK;
}

// Analysis_Average::Analyze()
Analysis::RetType Analysis_Average::Analyze() {
  if (writeHeader_)
    outfile_->Printf("%-6s %10s %10s %10s %10s %10s %10s %s\n",
                     "#Set", "Average", "Stdev", "Ymin", "YminIdx", 
                     "Ymax", "YmaxIdx", "Name");
  for (Array1D::const_iterator DS = input_dsets_.begin();
                               DS != input_dsets_.end(); ++DS)
  {
    if ( (*DS)->Size() < 1)
      mprintf("Warning: Set \"%s\" has no data.\n", (*DS)->legend());
    else {
      double Ymin = (*DS)->Dval(0);
      unsigned int idxYmin = 0;
      double Ymax = (*DS)->Dval(0);
      unsigned int idxYmax = 0;
      // TODO X min max? 
      double stdev = 0.0;
      double avg = (*DS)->Avg( stdev );
      // Find min/max and indices
      for (unsigned int idx = 1; idx != (*DS)->Size(); idx++) {
        double Yval = (*DS)->Dval(idx);
        if (Yval < Ymin) {
          Ymin = Yval;
          idxYmin = idx;
        }
        if (Yval > Ymax) {
          Ymax = Yval;
          idxYmax = idx;
        }
      }
      outfile_->Printf("%-6u %10.4g %10.4g %10.4g %10u %10.4g %10u \"%s\"\n",
                       DS - input_dsets_.begin(), avg, stdev, 
                       Ymin, idxYmin+1, Ymax, idxYmax+1, (*DS)->legend());
    }
  }
  return Analysis::OK;
}

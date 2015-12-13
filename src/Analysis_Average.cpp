#include <cmath>
#include "Analysis_Average.h"
#include "CpptrajStdio.h"

void Analysis_Average::Help() const {
  mprintf("\t<dset0> [<dset1> ...] [out <file>] [noheader] [oversets [name <outset>]]\n"
          "  Calculate the average, standard deviation, min, and max of given data sets.\n"
          "  If 'oversets' is specified calculate the average over all sets.\n");
}

// Analysis_Average::Setup()
Analysis::RetType Analysis_Average::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  calcAvgOverSets_ = analyzeArgs.hasKey("oversets");
  writeHeader_ = !analyzeArgs.hasKey("noheader");
  DataFile* setfile = 0;
  if (calcAvgOverSets_) {
    setfile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
    avgOfSets_ = setup.DSL().AddSet(DataSet::DOUBLE, analyzeArgs.GetStringKey("name"), "AVERAGE");
    if (avgOfSets_ == 0) return Analysis::ERR;
    sdOfSets_ = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(avgOfSets_->Meta().Name(), "SD"));
    if (sdOfSets_ == 0) return Analysis::ERR;
    if (setfile != 0) {
      setfile->AddDataSet( avgOfSets_ );
      setfile->AddDataSet( sdOfSets_ );
    }
  } else {
    outfile_ = setup.DFL().AddCpptrajFile(analyzeArgs.GetStringKey("out"), "DataSet Average",
                                     DataFileList::TEXT, true);
    if (outfile_ == 0) return Analysis::ERR;
  }
  // Select datasets from remaining args
  if (input_dsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), setup.DSL() )) {
    mprinterr("Error: Could not add data sets.\n");
    return Analysis::ERR;
  }
  if (input_dsets_.empty()) {
    mprinterr("Error: No input data sets.\n");
    return Analysis::ERR;
  }

  mprintf("    AVERAGE:");
  if (calcAvgOverSets_) {
    mprintf(" Calculating average over %i data sets.\n", input_dsets_.size());
    mprintf("\tData set base name '%s'", avgOfSets_->Meta().Name().c_str());
    if (setfile != 0) mprintf(", written to %s", setfile->DataFilename().full());
    mprintf("\n");
  } else {
    mprintf(" Calculating average of %i data sets.\n", input_dsets_.size());
    mprintf("\tWriting results to %s\n", outfile_->Filename().full());
  }
  //for (Array1D::const_iterator set = input_dsets_.begin(); set != input_dsets_.end(); ++set)
  //  mprintf("\t%s\n", (*set)->legend());

  return Analysis::OK;
}

// Analysis_Average::Analyze()
Analysis::RetType Analysis_Average::Analyze() {
  if (calcAvgOverSets_) {
    mprintf("\tCalculating average over sets:\n");
    // Ensure data sets are the same size.
    size_t Ndata = 0;
    for (Array1D::const_iterator DS = input_dsets_.begin(); DS != input_dsets_.end(); ++DS) {
      mprintf("\t%s\n", (*DS)->legend());
      if (DS == input_dsets_.begin())
        Ndata = (*DS)->Size();
      else if (Ndata != (*DS)->Size()) {
        mprinterr("Error: Set %s size %zu does not match first set size %zu\n",
                  (*DS)->legend(), (*DS)->Size(), Ndata);
        return Analysis::ERR;
      }
    }
    // Loop over all data
    double nsets = (double)input_dsets_.size();
    for (unsigned int i = 0; i != Ndata; i++) {
      double avg = 0.0;
      double sd = 0.0;
      for (Array1D::const_iterator DS = input_dsets_.begin(); DS != input_dsets_.end(); ++DS) {
        double dval = (*DS)->Dval( i );
        avg += dval;
        sd += (dval * dval);
      }
      avg /= nsets;
      sd /= nsets;
      sd -= (avg * avg);
      if (sd > 0.0)
        sd = sqrt( sd );
      else
        sd = 0.0;
      avgOfSets_->Add(i, &avg);
      sdOfSets_->Add(i, &sd);
    }
  } else {
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
  }
  return Analysis::OK;
}

#include <cmath>
#include "Analysis_Average.h"
#include "CpptrajStdio.h"

Analysis_Average::Analysis_Average() :
  avgOfSets_(0),
  sdOfSets_(0),
  data_avg_(0),
  data_sd_(0),
  data_ymin_(0),
  data_ymax_(0),
  data_yminIdx_(0),
  data_ymaxIdx_(0),
  data_names_(0),
  calcAvgOverSets_(false),
  toStdout_(false)
{}

void Analysis_Average::Help() const {
  mprintf("\t<dset0> [<dset1> ...] [torsion] [out <file>] [oversets]\n"
          "\t[name <output setname>] [nostdout]\n"
          "  Calculate the average, standard deviation, min, and max of given data sets.\n"
          "  If 'oversets' is specified calculate the average over all sets.\n"
          "  If 'nostdout' specified do not write averages to STDOUT when 'out' not specified.\n");
}

// Analysis_Average::Setup()
Analysis::RetType Analysis_Average::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  calcAvgOverSets_ = analyzeArgs.hasKey("oversets");
  DataFile* outfile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  toStdout_ = (outfile == 0 && !calcAvgOverSets_);
  if (analyzeArgs.hasKey("nostdout")) toStdout_ = false;
  bool convertToTorsion = analyzeArgs.hasKey("torsion");
  std::string dsname = analyzeArgs.GetStringKey("name");
  if (dsname.empty())
    dsname = setup.DSL().GenerateDefaultName("AVERAGE");
  // Select datasets from remaining args
  if (input_dsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), setup.DSL() )) {
    mprinterr("Error: Could not add data sets.\n");
    return Analysis::ERR;
  }
  if (input_dsets_.empty()) {
    mprinterr("Error: No input data sets.\n");
    return Analysis::ERR;
  }
  if (convertToTorsion) {
    for (Array1D::const_iterator set = input_dsets_.begin(); set != input_dsets_.end(); ++set) {
      MetaData md = (*set)->Meta();
      md.SetScalarMode( MetaData::M_TORSION );
      (*set)->SetMeta( md );
    }
  }
  // Set up output data sets for the specified calculation mode
  if (calcAvgOverSets_) {
    if (convertToTorsion)
      mprintf("Warning: 'torsion' keyword not used with 'oversets'\n");
    avgOfSets_ = setup.DSL().AddSet(DataSet::DOUBLE, dsname);
    if (avgOfSets_ == 0) return Analysis::ERR;
    sdOfSets_ = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "SD"));
    if (sdOfSets_ == 0) return Analysis::ERR;
    if (outfile != 0) {
      outfile->AddDataSet( avgOfSets_ );
      outfile->AddDataSet( sdOfSets_ );
    }
  } else {
    MetaData md(dsname, "avg");
    data_avg_ = setup.DSL().AddSet(DataSet::DOUBLE, md);
    md.SetAspect("sd");
    data_sd_ = setup.DSL().AddSet(DataSet::DOUBLE, md);
    md.SetAspect("ymin");
    data_ymin_ = setup.DSL().AddSet(DataSet::DOUBLE, md);
    md.SetAspect("ymax");
    data_ymax_ = setup.DSL().AddSet(DataSet::DOUBLE, md);
    md.SetAspect("yminidx");
    data_yminIdx_ = setup.DSL().AddSet(DataSet::INTEGER, md);
    md.SetAspect("ymaxidx");
    data_ymaxIdx_ = setup.DSL().AddSet(DataSet::INTEGER, md);
    md.SetAspect("names");
    data_names_ = setup.DSL().AddSet(DataSet::STRING, md);
    if (data_avg_ == 0 || data_sd_ == 0 || data_ymin_ == 0 || data_ymax_ == 0 ||
        data_yminIdx_ == 0 || data_ymaxIdx_ == 0 || data_names_ == 0)
      return Analysis::ERR;
    if (outfile != 0) {
      outfile->AddDataSet(data_avg_);
      outfile->AddDataSet(data_sd_);
      outfile->AddDataSet(data_ymin_);
      outfile->AddDataSet(data_ymax_);
      outfile->AddDataSet(data_yminIdx_);
      outfile->AddDataSet(data_ymaxIdx_);
      outfile->AddDataSet(data_names_);
    }
  }

  mprintf("    AVERAGE:");
  if (calcAvgOverSets_) {
    mprintf(" Calculating average over %zu data sets.\n", input_dsets_.size());
    mprintf("\tAverage stored in data set '%s'\n", avgOfSets_->legend());
    mprintf("\tStandard deviation stored in data set '%s'\n", sdOfSets_->legend());
  } else {
    mprintf(" Calculating average of %zu data sets.\n", input_dsets_.size());
    mprintf("\tData set base name '%s'\n", data_avg_->Meta().Name().c_str());
    if (convertToTorsion)
      mprintf("\tData sets will be marked as torsions\n");
  }
  if (outfile != 0) mprintf("\tOutput to to '%s'\n", outfile->DataFilename().full());
  //for (Array1D::const_iterator set = input_dsets_.begin(); set != input_dsets_.end(); ++set)
  //  mprintf("\t%s\n", (*set)->legend());

  return Analysis::OK;
}

// Analysis_Average::Analyze()
Analysis::RetType Analysis_Average::Analyze() {
  if (calcAvgOverSets_) {
    mprintf("\tCalculating average over sets:");
    // Ensure data sets are the same size.
    size_t Ndata = 0;
    for (Array1D::const_iterator DS = input_dsets_.begin(); DS != input_dsets_.end(); ++DS) {
      mprintf(" %s", (*DS)->legend());
      if (DS == input_dsets_.begin())
        Ndata = (*DS)->Size();
      else if (Ndata != (*DS)->Size()) {
        mprinterr("Error: Set %s size %zu does not match first set size %zu\n",
                  (*DS)->legend(), (*DS)->Size(), Ndata);
        return Analysis::ERR;
      }
    }
    mprintf("\n");
    // For now, set output set dimensions to first set dimension.
    avgOfSets_->SetDim(Dimension::X, input_dsets_[0]->Dim(0));
    sdOfSets_->SetDim(Dimension::X, input_dsets_[0]->Dim(0));
    // Loop over all data TODO average X as well?
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
    Dimension Xdim(1, 1, "Set");
    data_avg_->SetDim(Dimension::X, Xdim);
    data_sd_->SetDim(Dimension::X, Xdim);
    data_ymin_->SetDim(Dimension::X, Xdim);
    data_ymax_->SetDim(Dimension::X, Xdim);
    data_yminIdx_->SetDim(Dimension::X, Xdim);
    data_ymaxIdx_->SetDim(Dimension::X, Xdim);
    data_names_->SetDim(Dimension::X, Xdim);
    // Default to better format for very large/small numbers
    TextFormat Fmt(TextFormat::GDOUBLE, 10, 4);
    data_avg_->SetupFormat() = Fmt;
    data_sd_->SetupFormat() = Fmt;
    data_ymin_->SetupFormat() = Fmt;
    data_ymax_->SetupFormat() = Fmt;
    Fmt = TextFormat(TextFormat::INTEGER, 10);
    data_yminIdx_->SetupFormat() = Fmt;
    data_ymaxIdx_->SetupFormat() = Fmt;
    int set = 0;
    for (Array1D::const_iterator DS = input_dsets_.begin();
                                 DS != input_dsets_.end(); ++DS, ++set)
    {
      if ( (*DS)->Size() < 1)
        mprintf("Warning: Set \"%s\" has no data.\n", (*DS)->legend());
      else {
        mprintf("\t%i: %s\n", set+1, (*DS)->legend());
        std::string legend_with_quotes("\"" + (*DS)->Meta().Legend() + "\"");
        data_names_->Add( set, legend_with_quotes.c_str() );
        double Ymin = (*DS)->Dval(0);
        int idxYmin = 0;
        double Ymax = (*DS)->Dval(0);
        int idxYmax = 0;
        // TODO X min max? 
        double stdev = 0.0;
        double avg = (*DS)->Avg( stdev );
        data_avg_->Add( set, &avg );
        data_sd_->Add( set, &stdev );
        // Find min/max and indices
        for (int idx = 1; idx != (int)(*DS)->Size(); idx++) {
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
        idxYmin++;
        idxYmax++;
        data_ymin_->Add( set, &Ymin );
        data_ymax_->Add( set, &Ymax );
        data_yminIdx_->Add( set, &idxYmin );
        data_ymaxIdx_->Add( set, &idxYmax );
      }
    }
  }
  if (toStdout_) {
    DataFile OUT;
    OUT.SetupStdout(0); // TODO debug?
    OUT.AddDataSet( data_avg_ );
    OUT.AddDataSet( data_sd_ );
    OUT.AddDataSet( data_ymin_ );
    OUT.AddDataSet( data_ymax_ );
    OUT.AddDataSet( data_yminIdx_ );
    OUT.AddDataSet( data_ymaxIdx_ );
    OUT.AddDataSet( data_names_ );
    OUT.WriteDataOut();
  }
  return Analysis::OK;
}

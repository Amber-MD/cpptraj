#include "Analysis_RunningAvg.h"
#include "CpptrajStdio.h"
#include "DataSet_Mesh.h"

// CONSTRUCTOR
Analysis_RunningAvg::Analysis_RunningAvg() : cumulative_(false), window_(5) {}

void Analysis_RunningAvg::Help() const {
  mprintf("\t<dset1> [<dset2> ...] [name <dsetname>] [out <filename>]\n"
          "\t[ [cumulative] | [window <window>] ]\n"
          "  Calculate running average of data in selected data set(s)\n");
}

// Analysis_RunningAvg::Setup()
Analysis::RetType Analysis_RunningAvg::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  DataFile* outfile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  std::string setname = analyzeArgs.GetStringKey("name");
  cumulative_ = analyzeArgs.hasKey("cumulative");
  window_ = analyzeArgs.getKeyDouble("window", 5);

  // The remaining arguments are the data sets to take running averages of
  if (dsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), setup.DSL() )) {
    mprinterr("Error: runningavg: Could not add data sets.\n");
    return Analysis::ERR;
  }

  // If setname is empty, generate a default name
  if (setname.empty())
    setname = setup.DSL().GenerateDefaultName( "runningavg" );
  // Setup output datasets. Use XY Mesh so X can be avgd as well.
  int idx = 0;
  for (Array1D::const_iterator DS = dsets_.begin(); DS != dsets_.end(); ++DS) {
    DataSet* dsout = setup.DSL().AddSet(DataSet::XYMESH, MetaData(setname, idx++));
    if (dsout == 0)
      return Analysis::ERR;
    dsout->SetLegend( "RunAvg(" + (*DS)->Meta().Legend() + ")" );
    outputData_.push_back( dsout );
    if (outfile != 0) outfile->AddDataSet( dsout );
  }

  if (cumulative_)
    mprintf("    RUNNINGAVG: Calculating the cumulative running average for %zu data sets:\n",
            dsets_.size());
  else
    mprintf("    RUNNINGAVG: Calculating the running average for %zu data sets"
            " with a %d-element window:\n", dsets_.size(), window_);
  for (Array1D::const_iterator set = dsets_.begin(); set != dsets_.end(); ++set)
      mprintf("\t%s\n", (*set)->legend());
  if ( outfile != 0 )
    mprintf("\tOutfile name: %s\n", outfile->DataFilename().base());

  return Analysis::OK;
}

// Analysis_RunningAvg::Analyze()
Analysis::RetType Analysis_RunningAvg::Analyze() {
  std::vector<DataSet*>::const_iterator dsout = outputData_.begin();
  for (Array1D::const_iterator DS = dsets_.begin(); DS != dsets_.end(); DS++, ++dsout)
  {
    DataSet_1D const& data = *(*DS);
    DataSet_Mesh& out = static_cast<DataSet_Mesh&>( *(*dsout) );
    if (data.Size() < 2)
      mprintf("Warning: Set '%s' size is less than 2. Skipping.\n", data.legend());
    else {
      // If input data set X dim does not have default min/step, set them.
      if (cumulative_) {
        // Cumulative running average.
        mprintf("\t\tCalculating Cumulative Running Average for set %s\n", data.legend());
        double running_sum = 0.0;
        for (unsigned int i = 0; i < data.Size(); i++) {
          running_sum += data.Dval(i);
          out.AddXY( data.Xcrd(i), running_sum / (double)(i + 1) );
        }
      } else {
        mprintf("\t\tCalculating Running Average for set %s\n", data.legend());
        // Straight running average. Calculate initial average over window.
        double dwindow = (double)window_;
        double sumx = 0.0;
        double sumy = 0.0;
        for (int i = 0; i < window_; i++) {
          sumx += data.Xcrd( i );
          sumy += data.Dval( i );
        }
        out.AddXY( sumx / dwindow, sumy / dwindow );
        for (int i = 1; i < ((int)data.Size() - window_ + 1); i++) {
          int nextwin = i + window_ - 1;
          int prevwin = i - 1;
          sumx = data.Xcrd(nextwin) - data.Xcrd(prevwin) + sumx;
          sumy = data.Dval(nextwin) - data.Dval(prevwin) + sumy;
          out.AddXY( sumx / dwindow, sumy / dwindow );
        }
      }
      // Set output data set X dimension label.
      out.ModifyDim(Dimension::X).SetLabel("X");
    }
  }

  return Analysis::OK;
}

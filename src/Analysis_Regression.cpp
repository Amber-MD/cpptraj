#include "Analysis_Regression.h"
#include "CpptrajStdio.h"
#include "DataSet_Mesh.h"

Analysis_Regression::Analysis_Regression() : statsout_(0), nx_(0) {}

void Analysis_Regression::Help() const {
  mprintf("\t<dset0> [<dset1> ...] [name <name>] [nx <nxvals>]\n"
          "\t[out <filename>] [statsout <filename>]\n"
          "  Calculate linear regression lines for given data sets. If 'nx' is\n"
          "  specified the output data sets will have <nxvals> X values ranging\n"
          "  from the set minimum to maximum, otherwise X values from the set(s)\n"
          "  being fit will be used.\n");
}

// Analysis_Regression::Setup()
Analysis::RetType Analysis_Regression::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  nx_ = analyzeArgs.getKeyInt("nx", -1);
  if (nx_ > -1 && nx_ < 2) {
    mprinterr("Error: 'nx' must be greater than 1 if specified.\n");
    return Analysis::ERR;
  }
  DataFile* outfile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  statsout_ = setup.DFL().AddCpptrajFile(analyzeArgs.GetStringKey("statsout"),
                                    "Linear regression stats", DataFileList::TEXT, true);
  if (statsout_ == 0) return Analysis::ERR;
  std::string setname = analyzeArgs.GetStringKey("name");
  // Select datasets from remaining args
  if (input_dsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), setup.DSL() )) {
    mprinterr("Error: Could not add data sets.\n");
    return Analysis::ERR;
  }
  if (input_dsets_.empty()) {
    mprinterr("Error: No input data sets.\n");
    return Analysis::ERR;
  }
  // Setup output data sets TODO slope and intercept data sets
  int idx = 0;
  if ( input_dsets_.size() == 1 )
    idx = -1; // Only one input set, no need to refer to it by index
  // If setname is empty generate a default name
  if (setname.empty())
    setname = setup.DSL().GenerateDefaultName( "LR" );
  DataSet::DataType dtype;
  if (nx_ > 1)
    dtype = DataSet::DOUBLE;
  else
    dtype = DataSet::XYMESH;
  for ( Array1D::const_iterator DS = input_dsets_.begin();
                                DS != input_dsets_.end(); ++DS)
  {
    DataSet* dsout = setup.DSL().AddSet( dtype, MetaData(setname, idx++) );
    if (dsout==0) return Analysis::ERR;
    dsout->SetLegend( "LR(" + (*DS)->Meta().Legend() + ")" );
    output_dsets_.push_back( (DataSet_1D*)dsout );
    if (outfile != 0) outfile->AddDataSet( dsout );
  }

  mprintf("    REGRESSION: Calculating linear regression of %i data sets.\n",
          input_dsets_.size());
  if (outfile != 0)
    mprintf("\tFit line output to %s\n", outfile->DataFilename().full());
  mprintf("\tFit statistics output to %s\n", statsout_->Filename().full());
  if (nx_ > 1)
    mprintf("\tUsing %i X values from input set min to max\n", nx_);
  else
    mprintf("\tUsing X values from input sets\n");

  return Analysis::OK;
}

// Analysis_Regression::Analyze()
Analysis::RetType Analysis_Regression::Analyze() {
  int nerr = 0;
  //outfile_.Printf("#SetNum\tAverage\tStdev\tMin\tMax\tName\n");
  Array1D::const_iterator dsout = output_dsets_.begin();
  for (Array1D::const_iterator DS = input_dsets_.begin();
                               DS != input_dsets_.end();
                             ++DS, ++dsout)
  {
    if ( (*DS)->Size() < 2)
      mprintf("Warning: Set \"%s\" does not have enough data for regression (%zu points).\n", 
              (*DS)->legend(), (*DS)->Size());
    else {
      double slope, intercept, correl;
      mprintf("  %zu: %s\n", DS - input_dsets_.begin(), (*DS)->legend());
      if (!statsout_->IsStream())
        statsout_->Printf("#Stats for %s\n", (*DS)->legend());
      int err = (*DS)->LinearRegression( slope, intercept, correl, statsout_ );
      nerr += err;
      if (err == 0) {
        // Calculate fitted function
        if (nx_ < 2) {
          DataSet_Mesh& outMesh = static_cast<DataSet_Mesh&>( *(*dsout) );
          for (unsigned int i = 0; i < (*DS)->Size(); i++) {
            double x = (*DS)->Xcrd( i );
            outMesh.AddXY( x, slope * x + intercept );
          }
        } else {
          // Get min and max
          double xmin = (*DS)->Xcrd( 0 );
          double xmax = xmin;
          for (unsigned int i = 1; i < (*DS)->Size(); i++) {
            double xval = (*DS)->Xcrd( i );
            xmin = std::min( xmin, xval );
            xmax = std::max( xmax, xval );
          }
          double xstep = (xmax - xmin) / (double)(nx_ - 1);
          double xval = xmin;
          for (int i = 0; i < nx_; i++) {
            double yval = slope * xval + intercept;
            (*dsout)->Add(i, &yval);
            xval += xstep;
          }
          (*dsout)->SetDim( Dimension::X, Dimension(xmin, xstep, "X") );
        }
      }
    }
  }
  if (nerr > 0) return Analysis::ERR;
  return Analysis::OK;
}

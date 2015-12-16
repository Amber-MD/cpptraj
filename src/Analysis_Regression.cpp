#include "Analysis_Regression.h"
#include "CpptrajStdio.h"
#include "DataSet_Mesh.h"

Analysis_Regression::Analysis_Regression() : statsout_(0) {}

void Analysis_Regression::Help() const {
  mprintf("\t<dset0> [<dset1> ...] [name <name>] [out <filename>] [statsout <filename>]\n"
          "  Calculate linear regression lines for given data sets.\n");
}

// Analysis_Regression::Setup()
Analysis::RetType Analysis_Regression::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
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
  for ( Array1D::const_iterator DS = input_dsets_.begin();
                                DS != input_dsets_.end(); ++DS)
  {
    DataSet* dsout = setup.DSL().AddSet( DataSet::XYMESH, MetaData(setname, idx++) );
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
        DataSet_Mesh& outMesh = static_cast<DataSet_Mesh&>( *(*dsout) );
        for (unsigned int i = 0; i < (*DS)->Size(); i++) {
          double x = (*DS)->Xcrd( i );
          outMesh.AddXY( x, slope * x + intercept );
        }
      }
    }
  }
  if (nerr > 0) return Analysis::ERR;
  return Analysis::OK;
}

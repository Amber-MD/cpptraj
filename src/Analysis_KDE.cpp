#include <cmath> // pow
#include "Analysis_KDE.h"
#include "CpptrajStdio.h"
#include "DataSet_double.h"
#include "Constants.h" // TWOPI

Analysis_KDE::Analysis_KDE() : 
  data_(0), bandwidth_(0.0), output_(0), Kernel_(&Analysis_KDE::GaussianKernel) {}

void Analysis_KDE::Help() {
  mprintf("\t<dataset> [bandwidth <bw>] [out <file>] [name <dsname>]\n");
  mprintf("\t[min <min>] [max <max] [step <step>] [bins <bins>]\n");
}

Analysis::RetType Analysis_KDE::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  Dimension Xdim;
  if (analyzeArgs.Contains("min"))
    Xdim.SetMin( analyzeArgs.getKeyDouble("min", 0.0) );
  if (analyzeArgs.Contains("max"))
    Xdim.SetMax( analyzeArgs.getKeyDouble("max", 0.0) );
  Xdim.SetStep( analyzeArgs.getKeyDouble("step", -1.0) );
  Xdim.SetBins( analyzeArgs.getKeyInt("bins", -1) );
  if (Xdim.Step() < 0.0 && Xdim.Bins() < 0) {
    mprinterr("Error: Must set either bins or step.\n");
    return Analysis::ERR;
  }
  std::string setname = analyzeArgs.GetStringKey("name");
  bandwidth_ = analyzeArgs.getKeyDouble("bandwidth", -1.0);
  std::string outfilename = analyzeArgs.GetStringKey("out");
  DataFile* outfile = 0;
  if (!outfilename.empty())
    outfile = DFLin->AddDataFile( outfilename, analyzeArgs );

  // Get data set
  data_ = datasetlist->GetDataSet( analyzeArgs.GetStringNext() );
  if (data_ == 0) {
    mprinterr("Error: No data set or invalid data set name specified\n");
    return Analysis::ERR;
  }
  if (data_->Ndim() != 1) {
    mprinterr("Error: Only 1D data sets supported.\n");
    return Analysis::ERR;
  }
  
  // Output data set
  output_ = datasetlist->AddSet(DataSet::DOUBLE, setname, "kde");
  output_->SetDim(Dimension::X, Xdim);
  if (outfile != 0) outfile->AddSet( output_ );

  mprintf("    KDE: Using gaussian KDE to histogram set \"%s\"\n", data_->Legend().c_str());
  if (bandwidth_ < 0.0)
    mprintf("\tBandwidth will be estimated.\n");
  else
    mprintf("\tBandwidth= %f\n", bandwidth_);
  return Analysis::OK;
}

const double Analysis_KDE::ONE_OVER_ROOT_TWOPI = 1.0 / sqrt( TWOPI );

double Analysis_KDE::GaussianKernel(double u) const {
  return ( ONE_OVER_ROOT_TWOPI * exp( -0.5 * u * u ) );
}

Analysis::RetType Analysis_KDE::Analyze() {
  Dimension& Xdim = output_->Dim(0);
  DataSet_1D const& In = static_cast<DataSet_1D const&>( *data_ );
  // Set output set dimensions from input set if necessary.
  if (!Xdim.MinIsSet())
    Xdim.SetMin( In.Min() );
  if (!Xdim.MaxIsSet())
    Xdim.SetMax( In.Max() );
  if (Xdim.CalcBinsOrStep()) return Analysis::ERR;
  Xdim.PrintDim();

  // Allocate output set
  DataSet_double& Out = static_cast<DataSet_double&>( *output_ );
  Out.Resize( Xdim.Bins() );

  // Estimate bandwidth from normal distribution approximation if necessary.
  if (bandwidth_ < 0.0) {
    double stdev;
    In.Avg( stdev );
    double N_to_1_over_5 = pow( (double)In.Size(), (-1.0/5.0) );
    bandwidth_ = 1.06 * stdev * N_to_1_over_5;
    mprintf("\tDetermined bandwidth from normal distribution approximation: %f\n", bandwidth_);
  }
 
  // Loop over input data
  double total = 0.0;
  for (unsigned int i = 0; i < In.Size(); i++) {
    double val = In.Dval(i);
    double increment = 1.0;
    total += increment;
    // Apply kernel across histogram
    for (unsigned int j = 0; j < Out.Size(); j++)
      Out[j] += (increment * (this->*Kernel_)( (Xdim.Coord(j) - val) / bandwidth_ ));
  }

  // Normalize
  for (unsigned int j = 0; j < Out.Size(); j++)
    Out[j] /= (total * bandwidth_);

  return Analysis::OK;
}

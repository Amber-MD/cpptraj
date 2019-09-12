#include "Analysis_Integrate.h"
#include "CpptrajStdio.h"
#include "DataSet_Mesh.h"

Analysis_Integrate::Analysis_Integrate() {}

void Analysis_Integrate::Help() const {
  mprintf("\t<dset0> [<dset1> ...] [out <outfile>] [name <outsetname>]\n"
          "  Integrate given data sets.\n");
}

Analysis::RetType Analysis_Integrate::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  std::string setname = analyzeArgs.GetStringKey("name");
  DataFile* outfile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  // Select datasets from remaining args
  if (input_dsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), setup.DSL() )) {
    mprinterr("Error: Could not add data sets.\n");
    return Analysis::ERR;
  }
  if (input_dsets_.empty()) {
    mprinterr("Error: No input data sets.\n");
    return Analysis::ERR;
  }

  // Set up output datasets
  if (outfile != 0) {
    for (Array1D::const_iterator dsIn = input_dsets_.begin();
                                 dsIn != input_dsets_.end(); ++dsIn)
    {
      DataSet* ds = setup.DSL().AddSet(DataSet::XYMESH, setname, "Int");
      if (ds == 0) return Analysis::ERR;
      ds->SetLegend( "Int(" + (*dsIn)->Meta().Legend() + ")" );
      outfile->AddDataSet( ds );
      output_dsets_.push_back( (DataSet_Mesh*)ds );
    }
  }
  
  mprintf("    INTEGRATE: Calculating integral for %zu data sets.\n",
          input_dsets_.size());
  if (outfile != 0) {
    if (!setname.empty())
      mprintf("\tOutput set name: %s\n", setname.c_str());
    mprintf("\tOutfile name: %s\n", outfile->DataFilename().base());
  }
  if (debugIn > 0) {
    for (Array1D::const_iterator set = input_dsets_.begin(); set != input_dsets_.end(); ++set)
      mprintf("\t%s\n", (*set)->legend());
  }
  return Analysis::OK;
}

Analysis::RetType Analysis_Integrate::Analyze() {
  double sum;
  for (unsigned int idx = 0; idx != input_dsets_.size(); idx++) {
    DataSet_1D const& inSet = *(input_dsets_[idx]);
    if (inSet.Size() < 1) {
      mprintf("Warning: Set '%s' has no data.\n", inSet.legend());
    } else {
      if (!output_dsets_.empty()) {
        DataSet_Mesh& outSet = static_cast<DataSet_Mesh&>( *(output_dsets_[idx]) );
        sum = inSet.Integrate( DataSet_1D::TRAPEZOID, outSet.SetMeshX(), outSet.SetMeshY() );
        outSet.SetDim(Dimension::X, inSet.Dim(0));
      } else
        sum = inSet.Integrate( DataSet_1D::TRAPEZOID );
      mprintf("\tIntegral of %s is %g\n", inSet.legend(), sum);
    }
  }
  return Analysis::OK;
}

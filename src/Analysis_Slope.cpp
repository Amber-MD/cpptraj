#include "Analysis_Slope.h"
#include "CpptrajStdio.h"
#include "DataSet_Mesh.h"

// Analysis_Slope::Help()
void Analysis_Slope::Help() const {
  mprintf("\t<dset0> [<dset1> ...] [out <outfile>] [name <outsetname>]\n"
          "\t[type {forward|backward|central}]\n"
          "  Calculate finite difference (default forward) of given data sets.\n");
}

const char* Analysis_Slope::dTypeStr_[] = {"forward", "backward", "central"};

// Analysis_Slope::Setup()
Analysis::RetType Analysis_Slope::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  std::string setname = analyzeArgs.GetStringKey("name");
  DataFile* outfile = setup.DFL().AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  std::string dtypestr = analyzeArgs.GetStringKey("type");
  if (!dtypestr.empty()) {
    if      (dtypestr == "forward" ) diffType_ = DataSet_1D::FORWARD;
    else if (dtypestr == "backward") diffType_ = DataSet_1D::BACKWARD;
    else if (dtypestr == "central" ) diffType_ = DataSet_1D::CENTRAL;
    else {
      mprinterr("Error: Unrecognized type: '%s'\n", dtypestr.c_str());
      return Analysis::ERR;
    }
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

  // Set up output datasets
  for (Array1D::const_iterator dsIn = input_dsets_.begin();
                               dsIn != input_dsets_.end(); ++dsIn)
  {
    DataSet* ds = setup.DSL().AddSet(DataSet::XYMESH, setname, "Diff");
    if (ds == 0) return Analysis::ERR;
    ds->SetLegend( "Diff(" + (*dsIn)->Meta().Legend() + ")" );
    if (outfile != 0) outfile->AddDataSet( ds );
    output_dsets_.push_back( (DataSet_Mesh*)ds );
  }
  
  mprintf("    SLOPE: Calculating %s finite difference for %zu data sets.\n", 
          dTypeStr_[diffType_], input_dsets_.size());
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

// Analysis_Slope::Analyze()
Analysis::RetType Analysis_Slope::Analyze() {
  for (unsigned int idx = 0; idx != input_dsets_.size(); idx++) {
    DataSet_1D const& inSet = *(input_dsets_[idx]);
    mprintf("\t%s\n", inSet.legend());
    if (inSet.Size() < 1) {
      mprintf("Warning: Set '%s' has no data.\n", inSet.legend());
    } else {
      DataSet_Mesh& outSet = static_cast<DataSet_Mesh&>( *(output_dsets_[idx]) );
      inSet.FiniteDifference( diffType_, outSet.SetMeshX(), outSet.SetMeshY() );
      outSet.SetDim(Dimension::X, inSet.Dim(0));
    }
  }
  return Analysis::OK;
}

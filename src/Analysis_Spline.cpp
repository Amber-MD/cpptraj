#include "Analysis_Spline.h"
#include "CpptrajStdio.h"

Analysis_Spline::Analysis_Spline() :
  outfile_(0),
  meshsize_(0),
  meshmin_(0.0),
  meshmax_(0.0),
  xmin_(0.0),
  xstep_(0.0)
{}

void Analysis_Spline::Help() {
  mprintf("\t<dset0> [<dset1> ...] [out <outfile>] meshsize <x>]\n"
          "\t[meshmin <mmin>] meshmax <mmax> [min <xmin>] [step <xstep>]\n"
          "\tSpline the given data sets.\n");
}

Analysis::RetType Analysis_Spline::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  std::string setname = analyzeArgs.GetStringKey("name");
  outfile_ = DFLin->AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  meshsize_ = analyzeArgs.getKeyInt("meshsize", 0);
  if (meshsize_ < 3) {
    mprinterr("Error: meshsize must be specified and > 2\n");
    return Analysis::ERR;
  }
  meshmin_ = analyzeArgs.getKeyDouble("meshmin", 0.0);
  meshmax_ = analyzeArgs.getKeyDouble("meshmax", -1.0);
  if (meshmax_ < meshmin_) {
    mprinterr("Error: meshmax must be specified and > meshmin\n");
    return Analysis::ERR;
  }
  xmin_ = analyzeArgs.getKeyDouble("min", 0.0);
  xstep_ = analyzeArgs.getKeyDouble("step", 1.0);
  // Select datasets from remaining args
  if (input_dsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), *datasetlist )) {
    mprinterr("Error: Could not add data sets.\n");
    return Analysis::ERR;
  }
  if (input_dsets_.empty()) {
    mprinterr("Error: No input data sets.\n");
    return Analysis::ERR;
  }

  // Set up output datasets
  if (outfile_ != 0) {
    for (Array1D::const_iterator dsIn = input_dsets_.begin();
                                 dsIn != input_dsets_.end(); ++dsIn)
    {
      DataSet* ds = datasetlist->AddSet(DataSet::XYMESH, setname, "Spline");
      if (ds == 0) return Analysis::ERR;
      ds->SetLegend( "Spline(" + (*dsIn)->Legend() + ")" );
      outfile_->AddSet( ds );
      output_dsets_.push_back( (DataSet_Mesh*)ds );
    }
    outfile_->Dim(Dimension::X).SetMin( meshmin_ );
    double meshstep = (meshmax_ - meshmin_) / (double)meshsize_;
    outfile_->Dim(Dimension::X).SetStep( meshstep );
  }

  mprintf("    SPLINE: Applying cubic splining to %u data sets, mesh size is %i.\n",
          input_dsets_.size(), meshsize_);
  mprintf("\tMesh min=%f, Mesh max=%f, data set xmin=%f, data set step=%f\n",
          meshmin_, meshmax_, xmin_, xstep_);
  if (outfile_ != 0) {
    if (!setname.empty())
      mprintf("\tOutput set name: %s\n", setname.c_str());
    mprintf("\tOutfile name: %s\n", outfile_->DataFilename().base());
  }
  //for (Array1D::const_iterator set = input_dsets_.begin(); set != input_dsets_.end(); ++set)
  //  mprintf("\t%s\n", (*set)->Legend().c_str());
  return Analysis::OK;
}

// TODO: Each data set needs dimension information.
Analysis::RetType Analysis_Spline::Analyze() {
  for (unsigned int idx = 0; idx < input_dsets_.size(); idx++) {
    output_dsets_[idx]->CalculateMeshX(meshsize_, meshmin_, meshmax_);
    // Calculate mesh Y values.
    output_dsets_[idx]->SetSplinedMesh( *input_dsets_[idx], xmin_, xstep_ );
    // DEBUG
    //for (unsigned int i = 0; i < output_dsets_[idx]->Size(); i++)
    //  mprintf("\t%12.4f %12.4g\n", output_dsets_[idx]->X(i), output_dsets_[idx]->Y(i));
  }
  return Analysis::OK;
}

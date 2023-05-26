#include "Analysis_CalcDiffusion.h"
#include "CpptrajStdio.h"

// Analysis_CalcDiffusion::Help()
void Analysis_CalcDiffusion::Help() const {
  mprintf("\t[crdset <coords set>] [maxlag <maxlag>] [<mask>] [time <dt>]\n"
          "\t[<name>] [out <file>]\n");
}

// Analysis_CalcDiffusion::Setup()
Analysis::RetType Analysis_CalcDiffusion::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringKey("crdset");
  TgtTraj_ = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
  if (TgtTraj_ == 0) {
    mprinterr("Error: Could not locate COORDS set corresponding to '%s'\n",
              setname.c_str());
    Help();
    return Analysis::ERR;
  }
  maxlag_ = analyzeArgs.getKeyInt("maxlag", -1);
  double time_ = analyzeArgs.getKeyDouble("time", 1.0);
  DataFile* outfile = setup.DFL().AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
  // Mask
  if (mask1_.SetMaskString( analyzeArgs.GetMaskNext() )) {
    mprinterr("Error: Could not set mask string.\n");
    return Analysis::ERR;
  }
  // Add DataSets
  std::string dsname_ = analyzeArgs.GetStringNext();
  if (dsname_.empty())
    dsname_ = setup.DSL().GenerateDefaultName("Diff");
  avg_x_ = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "X"));
  avg_y_ = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "Y"));
  avg_z_ = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "Z"));
  avg_r_ = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "R"));
  avg_a_ = setup.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname_, "A"));
  if (avg_x_ == 0 || avg_y_ == 0 || avg_z_ == 0 || avg_r_ == 0 || avg_a_ == 0) {
    mprinterr("Error: Could not allocate one or more average diffusion sets.\n");
    return Analysis::ERR;
  }
  if (outfile != 0) {
    outfile->AddDataSet( avg_r_ );
    outfile->AddDataSet( avg_x_ );
    outfile->AddDataSet( avg_y_ );
    outfile->AddDataSet( avg_z_ );
    outfile->AddDataSet( avg_a_ );
  }
  // Set X dim
  Dimension Xdim_ = Dimension(0.0, time_, "Time");
  avg_x_->SetDim(Dimension::X, Xdim_);
  avg_y_->SetDim(Dimension::X, Xdim_);
  avg_z_->SetDim(Dimension::X, Xdim_);
  avg_r_->SetDim(Dimension::X, Xdim_);
  avg_a_->SetDim(Dimension::X, Xdim_);
  // Set up diffusion sets
  if (results_.CreateDiffusionSets(setup.DSL(), dsname_))
    return Analysis::ERR;

  mprintf("    CALCDIFFUSION: Calculating diffusion from COORDS set '%s'\n", TgtTraj_->legend());
  if (maxlag_ > 0)
    mprintf("\tMaximum lag is %i frames.\n", maxlag_);
  mprintf("\tUsing atoms selected by mask '%s'\n", mask1_.MaskString());

  return Analysis::OK;
}

// Analysis_CalcDiffusion::Analyze()
Analysis::RetType Analysis_CalcDiffusion::Analyze() {
  if (TgtTraj_->Size() < 1) {
    mprinterr("Error: COORDS set '%s' is empty.\n", TgtTraj_->legend());
    return Analysis::ERR;
  }
  if (TgtTraj_->Top().SetupIntegerMask( mask1_ )) {
    mprinterr("Error: Could not set up mask '%s'\n", mask1_.MaskString());
    return Analysis::ERR;
  }
  mask1_.MaskInfo();
  if (mask1_.None()) {
    mprinterr("Error: Nothing selected by mask '%s'\n", mask1_.MaskString());
    return Analysis::ERR;
  }

  return Analysis::OK;
}

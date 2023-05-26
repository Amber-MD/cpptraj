#include "Analysis_CalcDiffusion.h"
#include "CpptrajStdio.h"
#include "DataSet_double.h"
#include <algorithm> // std::min
#include <cmath> // sqrt

/** CONSTRUCTOR */
Analysis_CalcDiffusion::Analysis_CalcDiffusion() :
  TgtTraj_(0),
  maxlag_(0),
  time_(0),
  avg_x_(0),
  avg_y_(0),
  avg_z_(0),
  avg_r_(0),
  avg_a_(0)
{}

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
  time_ = analyzeArgs.getKeyDouble("time", 1.0);
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
  if (TgtTraj_->Size() == 1) {
    mprinterr("Error: COORDS set '%s' has only 1 frame.\n");
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

  int maxframes = (int)TgtTraj_->Size();
  if (maxlag_ < 1) {
    maxlag_ = (maxframes / 2);
  } else {
    int halfmax = (maxframes / 2);
    if (maxlag_ > halfmax)
      mprintf("Warning: Lag %i is more than half the number of input frames %i\n", maxlag_, halfmax);
  }
  int stopframe = maxframes - maxlag_;

  mprintf("\tCalculating diffusion from set '%s' for atoms in mask '%s' from t=0 to %g ps.\n",
          TgtTraj_->legend(), mask1_.MaskString(), (double)stopframe * time_);
  mprintf("\tMax lag is %i frames.\n", maxlag_);
  mprintf("\tUsing frames 1 to %i as time origins.\n", stopframe+1);

  if (stopframe < 1) {
    mprinterr("Error: Stop frame is less than 1.\n");
    return Analysis::ERR;
  }

  // Allocate sets
  DataSet_double& AX = static_cast<DataSet_double&>( *avg_x_ );
  DataSet_double& AY = static_cast<DataSet_double&>( *avg_y_ );
  DataSet_double& AZ = static_cast<DataSet_double&>( *avg_z_ );
  DataSet_double& AA = static_cast<DataSet_double&>( *avg_a_ );
  DataSet_double& AR = static_cast<DataSet_double&>( *avg_r_ );
  AX.Resize( maxlag_ );
  AY.Resize( maxlag_ );
  AZ.Resize( maxlag_ );
  AA.Resize( maxlag_ );
  AR.Resize( maxlag_ );

  std::vector<unsigned int> count( maxlag_, 0 );

  int idx0, idx1;
  Frame frm0 = TgtTraj_->AllocateFrame();
  Frame frm1 = frm0;

  // LOOP OVER FRAMES
  for (idx0 = 0; idx0 <= stopframe; idx0++)
  {
//    mprintf("DEBUG: (t=%g) %i to", (double)idx0*time_, idx0);
    TgtTraj_->GetFrame(idx0, frm0);
    int endidx = std::min(idx0 + maxlag_, maxframes);
    int tidx = 0;
    for (idx1 = idx0; idx1 < endidx; idx1++, tidx++)
    {
//      mprintf(" %i", idx1);
      // TODO for idx1==idx0 this is the same frame
      TgtTraj_->GetFrame(idx1, frm1);
      // Loop over atoms
      for (AtomMask::const_iterator at = mask1_.begin(); at != mask1_.end(); ++at)
      {
        const double* xyz0 = frm0.XYZ( *at );
        const double* xyz1 = frm1.XYZ( *at );
        double delx = xyz1[0] - xyz0[0];
        double dely = xyz1[1] - xyz0[1];
        double delz = xyz1[2] - xyz0[2];
        // Calc distances for this atom
        double distx = delx * delx;
        double disty = dely * dely;
        double distz = delz * delz;
        double dist2 = distx + disty + distz;
//        mprintf("DEBUG: At=%i  frm %i to %i  t=%g  d2=%g\n", *at+1, idx0+1, idx1+1, (double)tidx*time_, dist2);
        // Accumulate distances
        AX[tidx] += distx;
        AY[tidx] += disty;
        AZ[tidx] += distz;
        AR[tidx] += dist2;
        AA[tidx] += sqrt(dist2);
        count[tidx]++;
      } // END loop over atoms
    } // END inner loop
//    mprintf("\n");
  } // END outer loop

  // Calculate averages
  for (idx0 = 0; idx0 < maxlag_; idx0++) {
    mprintf("DEBUG: Average at t=%g is from %u data points.\n", (double)idx0*time_, count[idx0]);
    if (count[idx0] > 0) {
      AX[idx0] /= (double)count[idx0];
      AY[idx0] /= (double)count[idx0];
      AZ[idx0] /= (double)count[idx0];
      AR[idx0] /= (double)count[idx0];
      AA[idx0] /= (double)count[idx0];
    }
  }

  return Analysis::OK;
}

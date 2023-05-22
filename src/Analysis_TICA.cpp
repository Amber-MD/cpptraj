#include "Analysis_TICA.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
Analysis_TICA::Analysis_TICA() :
  TgtTraj_(0),
  lag_(0)
{}

// Analysis_TICA::Help()
void Analysis_TICA::Help() const {
  mprintf("[crdset <set name>] [lag <time lag>] [mask <mask>]\n");
}

// Analysis_TICA::Setup()
Analysis::RetType Analysis_TICA::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringKey("crdset");
  TgtTraj_ = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
  if (TgtTraj_ == 0) {
    mprinterr("Error: Could not locate COORDS set corresponding to %s\n",
              setname.c_str());
    Help();
    return Analysis::ERR;
  }
  // Other keywords
  lag_ = analyzeArgs.getKeyInt("lag", 1);
  std::string maskstr = analyzeArgs.GetStringKey("mask");
  if (mask1_.SetMaskString( maskstr )) {
    mprinterr("Error: Could not set atom mask string '%s'\n", maskstr.c_str());
    return Analysis::ERR;
  }

  // Print analysis info
  mprintf("    TICA: Time independent correlation analysis.\n");
  mprintf("\tUsing coordinates from set '%s'\n", TgtTraj_->legend());
  mprintf("\tUsing atoms selected by mask '%s'\n", mask1_.MaskString());
  mprintf("\tTime lag: %i frames.\n", lag_);

  return Analysis::OK;
}

// Analysis_TICA::Analyze()
Analysis::RetType Analysis_TICA::Analyze() {
  unsigned int Nframes = TgtTraj_->Size();
  if (Nframes < 1) {
    mprinterr("Error: No frames to analyze.\n");
    return Analysis::ERR;
  }
  // Evaluate mask
  if ( TgtTraj_->Top().SetupIntegerMask( mask1_ ) ) {
    mprinterr("Error: Could not evaluate atom mask '%s'\n", mask1_.MaskString());
    return Analysis::ERR;
  }
  mask1_.MaskInfo();
  if (mask1_.None()) {
    mprinterr("Error: No atoms selected by mask '%s'\n", mask1_.MaskString());
    return Analysis::ERR;
  }
  // Allocate frames
  Frame coords0;
  coords0.SetupFrameFromMask( mask1_, TgtTraj_->Top().Atoms(), TgtTraj_->CoordsInfo() );
  Frame coords1 = coords0;
  // Diagonal vectors
  typedef std::vector<Vec3> Varray;
  Varray vect(mask1_.Nselected(), Vec3(0.0));
  Varray vect2(mask1_.Nselected(), Vec3(0.0));
  // Loop over frames
  for (unsigned int frm0 = 0; frm0 < Nframes; frm0++) {
    mprintf("DEBUG: Frame %i\n", frm0);
    TgtTraj_->GetFrame(frm0, coords0, mask1_);
    // Covariance
    for (int idx2 = 0; idx2 < mask1_.Nselected(); idx2++) {
      Vec3 XYZj( coords0.XYZ(mask1_[idx2]) );
      vect[idx2] += XYZj;
      vect2[idx2] += XYZj.Squared();
    } // END outer loop over idx2
  }
  return Analysis::OK;
}

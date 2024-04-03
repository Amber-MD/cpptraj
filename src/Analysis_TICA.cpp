#include "Analysis_TICA.h"
#include "CoordCovarMatrix_Full.h"
#include "CoordCovarMatrix_Half.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"

/** CONSTRUCTOR */
Analysis_TICA::Analysis_TICA() :
  TgtTraj_(0),
  lag_(0),
  useMass_(false),
  debugC0_(0),
  debugCT_(0)
{
  SetHidden(true);
}

// Analysis_TICA::Help()
void Analysis_TICA::Help() const {
  mprintf("[crdset <set name>] [lag <time lag>] [mask <mask>] [mass]\n");
}

// Analysis_TICA::Setup()
Analysis::RetType Analysis_TICA::Setup(ArgList& analyzeArgs, AnalysisSetup& setup, int debugIn)
{
  // Attempt to get coords dataset or datasets from datasetlist
  TgtTraj_ = 0;
  sets_.clear();
  std::string setname = analyzeArgs.GetStringKey("crdset");
  std::string dataarg = analyzeArgs.GetStringKey("data");
  if (!setname.empty() && !dataarg.empty()) {
    mprinterr("Error: Specify either 'crdset' or 'data', not both.\n");
    return Analysis::ERR;
  }
  if (!setname.empty()) {
    TgtTraj_ = (DataSet_Coords*)setup.DSL().FindCoordsSet( setname );
    if (TgtTraj_ == 0) {
      mprinterr("Error: Could not locate COORDS set corresponding to %s\n",
                setname.c_str());
      Help();
      return Analysis::ERR;
    }
  } else if (!dataarg.empty()) {
    while (!dataarg.empty()) {
      if (sets_.AddSetsFromArgs( ArgList(dataarg), setup.DSL() )) {
        mprinterr("Error: Could not add data sets using argument '%s'\n", dataarg.c_str());
        return Analysis::ERR;
      }
      dataarg = analyzeArgs.GetStringKey("data");
    }
  } else {
    mprinterr("Error: Must specify either 'crdset' or 'data'.\n");
    return Analysis::ERR;
  }
  // Other keywords
  lag_ = analyzeArgs.getKeyInt("lag", 1);
  std::string maskstr = analyzeArgs.GetStringKey("mask");
  if (mask1_.SetMaskString( maskstr )) {
    mprinterr("Error: Could not set atom mask string '%s'\n", maskstr.c_str());
    return Analysis::ERR;
  }
  maskstr = analyzeArgs.GetStringKey("mask2");
  if (!maskstr.empty()) {
    mprintf("DEBUG: Second mask detected.\n");
    if (mask2_.SetMaskString( maskstr )) {
      mprinterr("Error: Could not set second atom mask string '%s'\n", maskstr.c_str());
      return Analysis::ERR;
    }
  }

  useMass_ = analyzeArgs.hasKey("mass");

  debugC0_ = setup.DFL().AddCpptrajFile(analyzeArgs.GetStringKey("debugc0"), "TICA C0 debug",
                                                                   DataFileList::TEXT, true);
  if (debugC0_ == 0) {
    mprinterr("Error: Could not open C0 debug file.\n");
    return Analysis::ERR;
  }
  debugCT_ = setup.DFL().AddCpptrajFile(analyzeArgs.GetStringKey("debugct"), "TICA CT debug",
                                                                   DataFileList::TEXT, true);
  if (debugCT_ == 0) {
    mprinterr("Error: Could not open CT debug file.\n");
    return Analysis::ERR;
  }

  // Print analysis info
  mprintf("    TICA: Time independent correlation analysis.\n");
  if (TgtTraj_ != 0) {
    mprintf("\tUsing coordinates from set '%s'\n", TgtTraj_->legend());
    mprintf("\tUsing atoms selected by mask '%s'\n", mask1_.MaskString());
    if (useMass_)
      mprintf("\tMass-weighted.\n");
    else
      mprintf("\tNot mass-weighted.\n");
  }
  if (!sets_.empty()) {
    mprintf("\tUsing %zu data sets:", sets_.size());
    for (Array1D::const_iterator it = sets_.begin(); it != sets_.end(); ++it)
      mprintf(" %s", (*it)->legend());
    mprintf("\n");
  }
  mprintf("\tTime lag: %i frames.\n", lag_);
  if (debugC0_ != 0)
    mprintf("\tDebug C0 output to %s\n", debugC0_->Filename().full());
  if (debugCT_ != 0)
    mprintf("\tDebug CT output to %s\n", debugCT_->Filename().full());

  return Analysis::OK;
}

// Analysis_TICA::Analyze()
Analysis::RetType Analysis_TICA::Analyze() {
  if (TgtTraj_ != 0)
    return analyze_crdset();
  else if (!sets_.empty())
    return analyze_datasets();

  return Analysis::ERR;
}

/** Analyze multiple 1D data sets. */
Analysis::RetType Analysis_TICA::analyze_datasets() {
  // Matrix - half
  CoordCovarMatrix_Half covarMatrix;
  if (covarMatrix.SetupMatrix( sets_.Array() )) {
    mprinterr("Error: Could not set up C0 matrix for data sets.\n");
    return Analysis::ERR;
  }
  covarMatrix.AddDataToMatrix( sets_.Array() );
  // Normalize
  if (covarMatrix.FinishMatrix()) {
    mprinterr("Error: Could not normalize coordinate covariance matrix for C0.\n");
    return Analysis::ERR;
  }
  // DEBUG PRINT
  covarMatrix.DebugPrint("C0", *debugC0_, " %12.6f");

  return Analysis::OK;
}

/** Analyze coordinates data set. */
Analysis::RetType Analysis_TICA::analyze_crdset() {
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
  // DEBUG
  if (mask2_.MaskStringSet()) {
    if ( TgtTraj_->Top().SetupIntegerMask( mask2_ ) ) {
      mprinterr("Error: Could not evaluate second atom mask '%s'\n", mask2_.MaskString());
      return Analysis::ERR;
    }
    mask2_.MaskInfo();
    if (mask2_.None()) {
      mprinterr("Error: No atoms selected by second mask '%s'\n", mask2_.MaskString());
      return Analysis::ERR;
    }
  }
  // Allocate frames
  Frame coords0 = TgtTraj_->AllocateFrame();
  //coords0.SetupFrameFromMask( mask1_, TgtTraj_->Top().Atoms(), TgtTraj_->CoordsInfo() );
  //Frame coords1 = coords0;
  // Matrix - half
  CoordCovarMatrix_Half covarMatrix;
  if (covarMatrix.SetupMatrix(TgtTraj_->Top().Atoms(), mask1_, useMass_ )) {
    mprinterr("Error: Could not set up C0 matrix for '%s'.\n", TgtTraj_->legend());
    return Analysis::ERR;
  }
  // Loop over frames
  for (unsigned int frm0 = 0; frm0 < Nframes; frm0++) {
    //mprintf("DEBUG: Frame %i\n", frm0);
    //TgtTraj_->GetFrame(frm0, coords0, mask1_);
    TgtTraj_->GetFrame(frm0, coords0);
    // Covariance
    covarMatrix.AddFrameToMatrix( coords0, mask1_ );
  } // END loop over frames

  // Normalize
  if (covarMatrix.FinishMatrix()) {
    mprinterr("Error: Could not normalize coordinate covariance matrix for C0.\n");
    return Analysis::ERR;
  }

  // DEBUG PRINT
  covarMatrix.DebugPrint("C0", *debugC0_);

  if (mask2_.MaskStringSet()) {
    // DEBUG
    // Matrix - full
    CoordCovarMatrix_Full CT;
    CT.SetupMatrix(TgtTraj_->Top().Atoms(), mask1_,
                   TgtTraj_->Top().Atoms(), mask2_, useMass_);
    // Loop over frames
    for (unsigned int frm0 = 0; frm0 < Nframes; frm0++) {
      //mprintf("DEBUG: Frame %i\n", frm0);
      //TgtTraj_->GetFrame(frm0, coords0, mask1_);
      TgtTraj_->GetFrame(frm0, coords0);
      // Covariance
      CT.AddFrameToMatrix( coords0, mask1_, coords0, mask2_ );
    } // END loop over frames

    // Normalize
    if (CT.FinishMatrix()) {
      mprinterr("Error: Could not normalize coordinate covariance matrix for CT.\n");
      return Analysis::ERR;
    }

    // DEBUG PRINT
    CT.DebugPrint("CT", *debugCT_);
  }

  return Analysis::OK;
}

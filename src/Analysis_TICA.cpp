#include "Analysis_TICA.h"
#include "CpptrajStdio.h"
#include "Matrix.h" // TODO DataSet?
#include <cmath> // sqrt

/** CONSTRUCTOR */
Analysis_TICA::Analysis_TICA() :
  TgtTraj_(0),
  lag_(0),
  useMass_(false)
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
  useMass_ = analyzeArgs.hasKey("mass");

  // Print analysis info
  mprintf("    TICA: Time independent correlation analysis.\n");
  mprintf("\tUsing coordinates from set '%s'\n", TgtTraj_->legend());
  mprintf("\tUsing atoms selected by mask '%s'\n", mask1_.MaskString());
  mprintf("\tTime lag: %i frames.\n", lag_);
  if (useMass_)
    mprintf("\tMass-weighted.\n");
  else
    mprintf("\tNot mass-weighted.\n");

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
  //Varray vect2(mask1_.Nselected(), Vec3(0.0));
  // Matrix - half
  Matrix<double> covarMatrix;
  covarMatrix.resize( mask1_.Nselected()*3, 0 );
  // Loop over frames
  for (unsigned int frm0 = 0; frm0 < Nframes; frm0++) {
    mprintf("DEBUG: Frame %i\n", frm0);
    TgtTraj_->GetFrame(frm0, coords0, mask1_);
    // Covariance
    Matrix<double>::iterator mat = covarMatrix.begin();
    for (int idx1 = 0; idx1 < mask1_.Nselected(); idx1++) {
      Vec3 XYZi( coords0.XYZ(idx1) );
      // Store veci and veci^2
      vect[idx1] += XYZi;
      //vect2[idx1] += XYZi.Squared();
      // Loop over X, Y, and Z of veci
      for (int iidx = 0; iidx < 3; iidx++) {
        double Vi = XYZi[iidx];
        // Diagonal
        for (int jidx = iidx; jidx < 3; jidx++)
          *(mat++) += Vi * XYZi[jidx]; // Vi * j{0,1,2}, Vi * j{1,2}, Vi * j{2}
        // Inner loop
        for (int idx2 = idx1 + 1; idx2 < mask1_.Nselected(); idx2++) {
          Vec3 XYZj( coords0.XYZ(idx2) );
          *(mat++) += Vi * XYZj[0];
          *(mat++) += Vi * XYZj[1];
          *(mat++) += Vi * XYZj[2];
        } // END inner loop over idx2
      } // END loop over x y z of veci
    } // END outer loop over idx1
  } // END loop over frames

  // Normalize
  double norm = 1.0 / (double)TgtTraj_->Size();
  for (Varray::iterator it = vect.begin(); it != vect.end(); ++it)
    *it *= norm;
  for (Matrix<double>::iterator it = covarMatrix.begin(); it != covarMatrix.end(); ++it)
    *it *= norm;
  // Calc <riri> - <ri><ri>
  //for (int k = 0; k < mask1_.Nselected(); k++) {
  //  vect2[k][0] -= (vect[k][0] * vect[k][0]);
  //  vect2[k][1] -= (vect[k][1] * vect[k][1]);
  //  vect2[k][2] -= (vect[k][2] * vect[k][2]);
  //}
  // Calc <rirj> - <ri><rj>
  double Mass = 1.0;
  double mass1 = 1.0;
  Matrix<double>::iterator mat = covarMatrix.begin();
  for (int idx1 = 0; idx1 < mask1_.Nselected(); idx1++) {
    if (useMass_)
      mass1 = coords0.Mass(idx1);
    for (int iidx = 0; iidx < 3; iidx++) {
      double Vi = vect[idx1][iidx];
      for (int idx2 = idx1; idx2 < mask1_.Nselected(); idx2++) {
        if (useMass_)
          Mass = sqrt( mass1 * coords0.Mass(idx2) );
        if (idx1 == idx2) {
          // Self
          for (int jidx = iidx; jidx < 3; jidx++) {
            *mat = (*mat - (Vi * vect[idx2][jidx])) * Mass;
            ++mat;
          }
        } else {
          for (int jidx = 0; jidx < 3; jidx++) {
            *mat = (*mat - (Vi * vect[idx2][jidx])) * Mass;
            ++mat;
          }
        }
      } // END inner loop over idx2
    } // END loop over elements of vect[idx1]
  } // END outer loop over idx1
  // DEBUG PRINT
  for (int row = 0; row < mask1_.Nselected()*3; row++) {
    for (int col = 0; col < mask1_.Nselected()*3; col++) {
      mprintf(" %6.3f", covarMatrix.element(col, row));
    }
    mprintf("\n");
  }
  return Analysis::OK;
}

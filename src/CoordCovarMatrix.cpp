#include "CoordCovarMatrix.h"
#include "AtomMask.h"
#include "Frame.h"

/** CONSTRUCTOR */
CoordCovarMatrix::CoordCovarMatrix() :
  nframes_(0)
{}

/** Clear the matrix */
void CoordCovarMatrix::Clear() {
  covarMatrix_.clear();
  vect_.clear();
  nframes_ = 0;
}

/** Add given Frame to the matrix. */
void CoordCovarMatrix::AddFrame(Frame const& frameIn, AtomMask const& maskIn)
{
  // Covariance
  MatType::iterator mat = covarMatrix_.begin();
  for (int idx1 = 0; idx1 < maskIn.Nselected(); idx1++) {
    Vec3 XYZi( frameIn.XYZ(idx1) );
    // Store veci and veci^2
    vect_[idx1] += XYZi;
    //vect2[idx1] += XYZi.Squared();
    // Loop over X, Y, and Z of veci
    for (int iidx = 0; iidx < 3; iidx++) {
      double Vi = XYZi[iidx];
      // Diagonal
      for (int jidx = iidx; jidx < 3; jidx++)
        *(mat++) += Vi * XYZi[jidx]; // Vi * j{0,1,2}, Vi * j{1,2}, Vi * j{2}
      // Inner loop
      for (int idx2 = idx1 + 1; idx2 < maskIn.Nselected(); idx2++) {
        Vec3 XYZj( frameIn.XYZ(idx2) );
        *(mat++) += Vi * XYZj[0];
        *(mat++) += Vi * XYZj[1];
        *(mat++) += Vi * XYZj[2];
      } // END inner loop over idx2
    } // END loop over x y z of veci
  } // END outer loop over idx1
  nframes_++;
}

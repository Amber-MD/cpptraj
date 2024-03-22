#include "CoordCovarMatrix.h"
#include "Atom.h"
#include "AtomMask.h"
#include "CpptrajStdio.h"
#include "Frame.h"

/** CONSTRUCTOR */
CoordCovarMatrix::CoordCovarMatrix() :
  nframes_(0),
  useMass_(false)
{}

/** Clear the matrix */
void CoordCovarMatrix::Clear() {
  covarMatrix_.clear();
  vect_.clear();
  mass_.clear();
  nframes_ = 0;
}

/** Set up array sizess and masses. */
int CoordCovarMatrix::SetupMatrix(std::vector<Atom> const& atoms,
                                  AtomMask const& maskIn, bool useMassIn)
{
  useMass_ = useMassIn;
  // TODO more size error checking
  nframes_ = 0;
  vect_.assign(maskIn.Nselected(), Vec3(0.0));
  //Varray vect2(mask1_.Nselected(), Vec3(0.0));
  // Matrix - half
  covarMatrix_.resize( maskIn.Nselected()*3, 0 );
  // Masses
  mass_.clear();
  if (useMassIn) {
    for (AtomMask::const_iterator at = maskIn.begin(); at != maskIn.end(); ++at)
      mass_.push_back( atoms[*at].Mass() );
  } else {
    for (int idx = 0; idx < maskIn.Nselected(); idx++)
      mass_.push_back( 1.0 );
  }
  return 0;
}

/** Add given Frame to the matrix. */
void CoordCovarMatrix::AddFrameToMatrix(Frame const& frameIn, AtomMask const& maskIn)
{
  // Covariance
  MatType::iterator mat = covarMatrix_.begin();
  for (int idx1 = 0; idx1 < maskIn.Nselected(); idx1++) {
    int at1 = maskIn[idx1];
    Vec3 XYZi( frameIn.XYZ(at1) );
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
        int at2 = maskIn[idx2];
        Vec3 XYZj( frameIn.XYZ(at2) );
        *(mat++) += Vi * XYZj[0];
        *(mat++) += Vi * XYZj[1];
        *(mat++) += Vi * XYZj[2];
      } // END inner loop over idx2
    } // END loop over x y z of veci
  } // END outer loop over idx1
  nframes_++;
}

/** Finish the matrix. */
int CoordCovarMatrix::FinishMatrix() {
  if (nframes_ < 1) {
    mprinterr("Error: No frames in coordinate covariance matrix.\n");
    return 1;
  }
  // Normalize
  double norm = 1.0 / (double)nframes_;
  for (Varray::iterator it = vect_.begin(); it != vect_.end(); ++it)
    *it *= norm;
  for (MatType::iterator it = covarMatrix_.begin(); it != covarMatrix_.end(); ++it)
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
  MatType::iterator mat = covarMatrix_.begin();
  for (unsigned int idx1 = 0; idx1 < mass_.size(); idx1++) {
    mass1 = mass_[idx1];
    for (int iidx = 0; iidx < 3; iidx++) {
      double Vi = vect_[idx1][iidx];
      for (unsigned int idx2 = idx1; idx2 < mass_.size(); idx2++) {
        Mass = sqrt( mass1 * mass_[idx2] );
        if (idx1 == idx2) {
          // Self
          for (int jidx = iidx; jidx < 3; jidx++) {
            *mat = (*mat - (Vi * vect_[idx2][jidx])) * Mass;
            ++mat;
          }
        } else {
          for (int jidx = 0; jidx < 3; jidx++) {
            *mat = (*mat - (Vi * vect_[idx2][jidx])) * Mass;
            ++mat;
          }
        }
      } // END inner loop over idx2
    } // END loop over elements of vect_[idx1]
  } // END outer loop over idx1
  return 0;
}

/** Debug print to STDOUT */
void CoordCovarMatrix::DebugPrint(const char* desc) const {
  if (desc != 0)
    mprintf("DEBUG: CoordCovarMatrix: %s\n", desc);
  for (unsigned int row = 0; row < covarMatrix_.Nrows(); row++) {
    for (unsigned int col = 0; col < covarMatrix_.Ncols(); col++) {
      mprintf(" %6.3f", covarMatrix_.element(col, row));
    }
    mprintf("\n");
  }
}

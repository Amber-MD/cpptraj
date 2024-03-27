#include "CoordCovarMatrix_Full.h"
#include "AtomMask.h"
#include "CpptrajStdio.h"
#include "Frame.h"

/** CONSTRUCTOR */
CoordCovarMatrix_Full::CoordCovarMatrix_Full()
{}

/** Clear the matrix */
void CoordCovarMatrix_Full::clearMat() {
//  vect2_1_.clear();
//  vect2_2_.clear();
  vect_1_.clear();
  vect_2_.clear();
  mass1_.clear();
  mass2_.clear();
}

/** Set up array sizess and masses. */
int CoordCovarMatrix_Full::SetupMatrix(std::vector<Atom> const& atoms1,
                                       AtomMask const& maskIn1,
                                       std::vector<Atom> const& atoms2,
                                       AtomMask const& maskIn2, bool useMassIn)
{
  // Matrix - full
  covarMatrix_.resize( maskIn1.Nselected()*3, maskIn2.Nselected()*3 );

  vect_1_.assign(maskIn1.Nselected(), Vec3(0.0));
//  vect2_1_.assign(maskIn1.Nselected(), Vec3(0.0));
  vect_2_.assign(maskIn2.Nselected(), Vec3(0.0));
//  vect2_2_.assign(maskIn2.Nselected(), Vec3(0.0));
  set_mass_array( mass1_, atoms1, maskIn1, useMassIn );
  set_mass_array( mass2_, atoms2, maskIn2, useMassIn );

  return 0;
}

/** Store diagonal (average and average^2) */
static inline void store_diagonal(std::vector<Vec3>& vect,
//                                  std::vector<Vec3>& vect2,
                                  Frame const& frameIn,
                                  AtomMask const& maskIn)
{
  for (int idx = 0; idx < maskIn.Nselected(); idx++)
  {
    const double* XYZ = frameIn.XYZ( maskIn[idx] );
    for (int ii = 0; ii < 3; ii++) {
      vect[idx][ii]  += XYZ[ii];
//      vect2[idx][ii] += (XYZ[ii] * XYZ[ii]);
    }
  }
}

/** Add selected atoms in given Frames to the matrix. */
void CoordCovarMatrix_Full::AddFrameToMatrix(Frame const& frameIn1, AtomMask const& maskIn1,
                                             Frame const& frameIn2, AtomMask const& maskIn2)
{
  //MatType::iterator mat = covarMatrix_.begin();
  for (int idx2 = 0; idx2 < maskIn2.Nselected(); idx2++)
  {
    Matrix<double>::iterator mat = covarMatrix_.begin() + ((idx2*3)*covarMatrix_.Ncols()); // NX
    const double* XYZj = frameIn2.XYZ( maskIn2[idx2] );
    for (int ny = 0; ny < 3; ny++) {
      double Vj = XYZj[ny];
      for (int idx1 = 0; idx1 < maskIn1.Nselected(); idx1++)
      {
        const double* XYZi  = frameIn1.XYZ( maskIn1[idx1] );
        *(mat++) += Vj * XYZi[0];
        *(mat++) += Vj * XYZi[1];
        *(mat++) += Vj * XYZi[2];
      }
    }
  }
  // Mask1/mask2 diagonal
//  store_diagonal(vect_1_, vect2_1_, frameIn1, maskIn1);
//  store_diagonal(vect_2_, vect2_2_, frameIn2, maskIn2);
  store_diagonal(vect_1_, frameIn1, maskIn1);
  store_diagonal(vect_2_, frameIn2, maskIn2);
  nframes_++;
}

/** Calculate <v^2> - <v><v> */
/*static inline void vect2_minus_vect(std::vector<Vec3>& vect2, std::vector<Vec3> const& vect)
{
  // Sanity check
  if (vect2.size() != vect.size()) {
    mprinterr("Internal Error: CoordCovarMatrix_Full: Vect2 size (%zu) != Vect size (%zu)\n",
              vect2.size(), vect.size());
    return;
  }
  for (unsigned int idx = 0; idx != vect2.size(); idx++) {
    Vec3& V2 = vect2[idx];
    Vec3 const& V1 = vect[idx];
    for (int ny = 0; ny < 3; ny++)
      V2[ny] -= (V1[ny] * V1[ny]);
  }
}*/

/** Finish processing covariance matrix */
int CoordCovarMatrix_Full::FinishMatrix() {
  if (nframes_ < 1) {
    mprinterr("Error: No frames in coordinate covariance matrix.\n");
    return 1;
  }
  // Normalize
  double norm = 1.0 / (double)nframes_;
  for (Varray::iterator it = vect_1_.begin(); it != vect_1_.end(); ++it)
    *it *= norm;
  for (Varray::iterator it = vect_2_.begin(); it != vect_2_.end(); ++it)
    *it *= norm;
//  for (Varray::iterator it = vect2_1_.begin(); it != vect2_1_.end(); ++it)
//    *it *= norm;
//  for (Varray::iterator it = vect2_2_.begin(); it != vect2_2_.end(); ++it)
//    *it *= norm;
  for (MatType::iterator it = covarMatrix_.begin(); it != covarMatrix_.end(); ++it)
    *it *= norm;
  mprintf("DEBUG: First 3 elements: %f %f %f\n",
          covarMatrix_.element(0, 0),
          covarMatrix_.element(1, 0),
          covarMatrix_.element(2, 0));
  mprintf("DEBUG: First 3 elements of V1: %f %f %f\n",
          vect_1_[0][0],
          vect_1_[0][1],
          vect_1_[0][2]);
  mprintf("DEBUG: First 3 elements of V2: %f %f %f\n",
          vect_2_[0][0],
          vect_2_[0][1],
          vect_2_[0][2]);
  // Calc <riri> - <ri><ri>
//  vect2_minus_vect(vect2_1_, vect_1_);
//  vect2_minus_vect(vect2_2_, vect_2_);
  // Calc <rirj> - <rj><rj>
  Matrix<double>::iterator mat = covarMatrix_.begin();
  for (unsigned int idx2 = 0; idx2 < vect_2_.size(); idx2++)
  {
    double mass2 = mass2_[idx2];
    Vec3 const& V2 = vect_2_[idx2];
    for (int ny = 0; ny < 3; ny++) {
      for (unsigned int idx1 = 0; idx1 < vect_1_.size(); idx1++) {
        double Mass = sqrt(mass2 * mass1_[idx1]);
        for (int nx = 0; nx < 3; nx++) {
          mprintf("mat = (%f - (%f * %f)) * %f\n", *mat, V2[ny], vect_1_[idx1][nx], Mass);
          *mat = (*mat - (V2[ny] * vect_1_[idx1][nx])) * Mass;
          ++mat;
        }
      } // END loop over vect_1_
    }
  } // END loop over vect_2_
  return 0;
}

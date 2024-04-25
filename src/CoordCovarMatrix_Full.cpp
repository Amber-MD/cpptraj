#include "CoordCovarMatrix_Full.h"
#include "AtomMask.h"
#include "CpptrajStdio.h"
#include <cmath> //sqrt

/** CONSTRUCTOR */
CoordCovarMatrix_Full::CoordCovarMatrix_Full()
{}

/** Clear the matrix */
void CoordCovarMatrix_Full::clearMat() {
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
  nelt_ = 3; // xyz
  unsigned int arraySize1 = (unsigned int)maskIn1.Nselected() * nelt_;
  unsigned int arraySize2 = (unsigned int)maskIn2.Nselected() * nelt_;
  // Matrix - full
  covarMatrix_.resize( arraySize1, arraySize2 );

  vect_1_.assign(arraySize1, 0);
  vect_2_.assign(arraySize2, 0);
  set_mass_array( mass1_, atoms1, maskIn1, useMassIn );
  set_mass_array( mass2_, atoms2, maskIn2, useMassIn );

  return 0;
}

/** Store diagonal (average and average^2) */
/*static inline void store_diagonal(std::vector<double>& vect,
//                                  std::vector<Vec3>& vect2,
                                  unsigned int nelt,
                                  Frame const& frameIn,
                                  AtomMask const& maskIn)
{
  for (int idx = 0; idx < maskIn.Nselected(); idx++)
  {
    const double* XYZ = frameIn.XYZ( maskIn[idx] );
    int eidx = idx * nelt;
    for (unsigned int ii = 0; ii < nelt; ii++) {
      vect[eidx+ii]  += XYZ[ii];
//      vect2[idx][ii] += (XYZ[ii] * XYZ[ii]);
    }
  }
}*/

/** Store average */
static inline void store_diagonal(std::vector<double>& vect, std::vector<double> const& arrayIn)
{
  for (unsigned int idx = 0; idx < arrayIn.size(); ++idx)
    vect[idx] += arrayIn[idx];
}

/** Add selected atoms in given Frames to the matrix. */
void CoordCovarMatrix_Full::AddFrameToMatrix(Frame const& frameIn1, AtomMask const& maskIn1,
                                             Frame const& frameIn2, AtomMask const& maskIn2)
{
  Darray array1, array2; // TODO make class vars
  get_frame_coords(array1, frameIn1, maskIn1);
  get_frame_coords(array2, frameIn2, maskIn2);
  AddToMatrix(array1, array2);
}

/** Add elements to the matrix. */
void CoordCovarMatrix_Full::AddToMatrix(Darray const& array1, Darray const& array2) {
  // sanity checks
  if (!has_valid_size(array1)) {
    mprinterr("Internal Error: CoordCovarMatrix_Full::AddToMatrix(): Incoming array1 size %zu not divisible by %u\n",
              array1.size(), nelt_);
    return;
  }
  if (!has_valid_size(array2)) {
    mprinterr("Internal Error: CoordCovarMatrix_Full::AddToMatrix(): Incoming array2 size %zu not divisible by %u\n",
              array2.size(), nelt_);
    return;
  }

  //MatType::iterator mat = covarMatrix_.begin();
  for (unsigned int idx2 = 0; idx2 < mass2_.size(); idx2++)
  {
    unsigned int eidx2 = idx2 * nelt_;
    Matrix<double>::iterator mat = covarMatrix_.begin() + (eidx2*covarMatrix_.Ncols()); // NX
    const double* XYZj = (&array2[0]) + eidx2;
    for (unsigned int ny = 0; ny < nelt_; ny++) {
      double Vj = XYZj[ny];
      for (unsigned int idx1 = 0; idx1 < mass1_.size(); idx1++)
      {
        unsigned int eidx1 = idx1 * nelt_;
        const double* XYZi = (&array1[0]) + eidx1;
        for (unsigned int nx = 0; nx < nelt_; nx++) {
          *(mat++) += Vj * XYZi[nx];
        }
      }
    }
  }
  // Mask1/mask2 diagonal
//  store_diagonal(vect_1_, vect2_1_, frameIn1, maskIn1);
//  store_diagonal(vect_2_, vect2_2_, frameIn2, maskIn2);
//  store_diagonal(vect_1_, nelt_, frameIn1, maskIn1);
//  store_diagonal(vect_2_, nelt_, frameIn2, maskIn2);
  store_diagonal(vect_1_, array1);
  store_diagonal(vect_2_, array2);
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
  for (Darray::iterator it = vect_1_.begin(); it != vect_1_.end(); ++it)
    *it *= norm;
  for (Darray::iterator it = vect_2_.begin(); it != vect_2_.end(); ++it)
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
          vect_1_[0],
          vect_1_[1],
          vect_1_[2]);
  mprintf("DEBUG: First 3 elements of V2: %f %f %f\n",
          vect_2_[0],
          vect_2_[1],
          vect_2_[2]);
  // Calc <riri> - <ri><ri>
//  vect2_minus_vect(vect2_1_, vect_1_);
//  vect2_minus_vect(vect2_2_, vect_2_);
  // Calc <rirj> - <rj><rj>
  Matrix<double>::iterator mat = covarMatrix_.begin();
  for (unsigned int idx2 = 0; idx2 < mass2_.size(); idx2++)
  {
    double mass2 = mass2_[idx2];
    unsigned int eidx2 = idx2 * nelt_;
    for (unsigned int ny = 0; ny < nelt_; ny++) {
      double V2 = vect_2_[eidx2+ny];
      for (unsigned int idx1 = 0; idx1 < mass1_.size(); idx1++) {
        double Mass = sqrt(mass2 * mass1_[idx1]);
        unsigned int eidx1 = idx1 * nelt_;
        for (unsigned int nx = 0; nx < nelt_; nx++) {
          double V1 = vect_1_[eidx1+nx];
          mprintf("mat = (%f - (%f * %f)) * %f\n", *mat, V2, V1, Mass);
          *mat = (*mat - (V2 * V1)) * Mass;
          ++mat;
        }
      } // END loop over vect_1_
    }
  } // END loop over vect_2_
  return 0;
}

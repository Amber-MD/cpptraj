#include "CoordCovarMatrix_Half.h"
#include "AtomMask.h"
#include "Constants.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"
#include <cmath> // sqrt

/** CONSTRUCTOR */
CoordCovarMatrix_Half::CoordCovarMatrix_Half()
{}

/** Clear the matrix */
void CoordCovarMatrix_Half::clearMat() {
  vect_.clear();
  mass_.clear();
}

/** Set up array for incoming data sets. */
int CoordCovarMatrix_Half::SetupMatrix(DSarray const& sets)
{
  int is_periodic = -1;
  
  for (DSarray::const_iterator it = sets.begin(); it != sets.end(); ++it)
  {
    if ((*it)->Meta().IsTorsionArray()) {
      if (is_periodic == -1)
        is_periodic = 1;
      else if (is_periodic != 1) {
        mprinterr("Error: If one set is not periodic all must be. Set '%s' is periodic.\n",
                  (*it)->legend());
        return 1;
      }
    } else {
      if (is_periodic == -1)
        is_periodic = 0;
      else if (is_periodic != 0) {
        mprinterr("Error: If one set is periodic all must be. Set '%'s is not periodic.\n",
                  (*it)->legend());
        return 1;
      }
    }
  }
  if (is_periodic) {
    mprintf("\tPeriodic data sets detected.\n");
    nelt_ = 2;
  } else {
    mprintf("\tNon-periodic data sets detected.\n");
    nelt_ = 1;
  }
  unsigned int matSize = sets.size() * nelt_;
  // Matrix - half
  covarMatrix_.resize( matSize, 0 );

  vect_.assign( matSize, 0 );

  mass_.resize( sets.size(), 1.0 );

  return 0;
}

/** Set up array sizess and masses. */
int CoordCovarMatrix_Half::SetupMatrix(std::vector<Atom> const& atoms,
                                       AtomMask const& maskIn, bool useMassIn)
{
  nelt_ = 3; // xyz
  unsigned int arraySize = (unsigned int)maskIn.Nselected() * nelt_;
  // Matrix - half
  covarMatrix_.resize( arraySize, 0 );

  vect_.assign( arraySize, 0 );

  set_mass_array(mass_, atoms, maskIn, useMassIn);

  return 0;
}

/** Add data from sets to the matrix. */
void CoordCovarMatrix_Half::AddDataToMatrix(DSarray const& sets)
{
  // TODO check empty input array
  // Check that sets have same size
  unsigned int maxFrames = sets.front()->Size();
  for (DSarray::const_iterator it = sets.begin(); it != sets.end(); ++it)
  {
    if ((*it)->Size() != maxFrames) {
      mprinterr("Error: Set '%s' does not have same size (%zu) as first set (%u)\n",
                (*it)->legend(), (*it)->Size(), maxFrames);
      return;
    }
  }
  Darray arrayIn;
  arrayIn.resize( sets.size() * nelt_ );
  if (nelt_ == 2) {
    for (unsigned int idx = 0; idx < maxFrames; idx++) {
      unsigned int jdx = 0;
      for (DSarray::const_iterator it = sets.begin(); it != sets.end(); ++it, jdx += 2)
      {
        double dval = (*it)->Dval(idx) * Constants::DEGRAD;
        arrayIn[jdx  ] = cos( dval );
        arrayIn[jdx+1] = sin( dval );
      }
      AddToMatrix(arrayIn);
    }
  } else if (nelt_ == 1) {
    for (unsigned int idx = 0; idx < maxFrames; idx++) {
      for (unsigned int jdx = 0; jdx < sets.size(); jdx++) {
        arrayIn[jdx] = sets[jdx]->Dval(idx);
      }
      AddToMatrix(arrayIn);
    }
  } else {
    // Sanity check
    mprinterr("Internal Error: CoordCovarMatrix_Half::AddDataToMatrix(): Unsupported nelt %u\n", nelt_);
  }
}

/** Add selected atoms in given Frame to the matrix. */
void CoordCovarMatrix_Half::AddFrameToMatrix(Frame const& frameIn, AtomMask const& maskIn)
{
  Darray arrayIn; // TODO class array?
  get_frame_coords(arrayIn, frameIn, maskIn);
  AddToMatrix(arrayIn);
}

/** Add elements to the matrix. */
void CoordCovarMatrix_Half::AddToMatrix(Darray const& arrayIn) {
  // sanity check
  if (!has_valid_size(arrayIn)) {
    mprinterr("Internal Error: CoordCovarMatrix_Half::AddToMatrix(): Incoming array size %zu not divisible by %u\n",
              arrayIn.size(), nelt_);
    return;
  }
  // Covariance
  MatType::iterator mat = covarMatrix_.begin();
  for (unsigned int idx2 = 0; idx2 < mass_.size(); idx2++) {
    unsigned int eidx2 = idx2 * nelt_;
    const double* XYZj = (&arrayIn[0]) + eidx2;
    // Store average 
    for (unsigned int ej = 0; ej < nelt_; ej++)
      vect_[eidx2+ej] += XYZj[ej];
    //XYZj.Print("XYZj");
    // Loop over X, Y, and Z of vecj
    for (unsigned int jidx = 0; jidx < nelt_; jidx++) {
      double Vj = XYZj[jidx];
      // Diagonal
      for (unsigned int ej = jidx; ej < nelt_; ej++)
        *(mat++) += Vj * XYZj[ej]; // Vj * j{0,1,2}, Vj * j{1,2}, Vj * j{2}
      // Inner loop
      for (unsigned int idx1 = idx2 + 1; idx1 < mass_.size(); idx1++) {
        unsigned int eidx1 = idx1 * nelt_;
        const double* XYZi = (&arrayIn[0]) + eidx1;
        for (unsigned int ei = 0; ei < nelt_; ei++) {
          *(mat++) += Vj * XYZi[ei];
        }
      } // END inner loop over idx1
    } // END loop over x y z of vecj
  } // END outer loop over idx2
  nframes_++;
}

/** Add given Frame to the matrix. */
/*void CoordCovarMatrix_Half::AddFrameToMatrix(Frame const& frameIn)
{
  // Covariance
  MatType::iterator mat = covarMatrix_.begin();
  for (int idx1 = 0; idx1 < frameIn.Natom(); idx1++) {
    Vec3 XYZi( frameIn.XYZ(idx1) );
    // Store veci and veci^2
    vect_[idx1] += XYZi;
    //XYZi.Print("XYZi");
    //vect2[idx1] += XYZi.Squared();
    // Loop over X, Y, and Z of veci
    for (int iidx = 0; iidx < 3; iidx++) {
      double Vi = XYZi[iidx];
      // Diagonal
      for (int jidx = iidx; jidx < 3; jidx++)
        *(mat++) += Vi * XYZi[jidx]; // Vi * j{0,1,2}, Vi * j{1,2}, Vi * j{2}
      // Inner loop
      for (int idx2 = idx1 + 1; idx2 < frameIn.Natom(); idx2++) {
        Vec3 XYZj( frameIn.XYZ(idx2) );
        *(mat++) += Vi * XYZj[0];
        *(mat++) += Vi * XYZj[1];
        *(mat++) += Vi * XYZj[2];
      } // END inner loop over idx2
    } // END loop over x y z of veci
  } // END outer loop over idx1
  nframes_++;
}*/

/** Finish the matrix. */
int CoordCovarMatrix_Half::FinishMatrix() {
  if (nframes_ < 1) {
    mprinterr("Error: No frames in coordinate covariance matrix.\n");
    return 1;
  }
  // Normalize
  double norm = 1.0 / (double)nframes_;
  for (Darray::iterator it = vect_.begin(); it != vect_.end(); ++it)
    *it *= norm;
  for (MatType::iterator it = covarMatrix_.begin(); it != covarMatrix_.end(); ++it)
    *it *= norm;
  // DEBUG print mean
  CpptrajFile outfile;
  outfile.OpenWrite("debug.mean.dat");
  for (Darray::iterator it = vect_.begin(); it != vect_.end(); ++it)
    outfile.Printf("%12.6f\n", *it);
  outfile.CloseFile();
  // Calc <riri> - <ri><ri>
  //for (int k = 0; k < mask1_.Nselected(); k++) {
  //  vect2[k][0] -= (vect[k][0] * vect[k][0]);
  //  vect2[k][1] -= (vect[k][1] * vect[k][1]);
  //  vect2[k][2] -= (vect[k][2] * vect[k][2]);
  //}
  // Calc <rirj> - <ri><rj>
  MatType::iterator mat = covarMatrix_.begin();
//  double TwoN = (double)( covarMatrix_->Ncols() * 2 );
  for (unsigned int idx2 = 0; idx2 < mass_.size(); idx2++) {
    double mass2 = mass_[idx2];
    for (unsigned int jidx = 0; jidx < nelt_; jidx++) {
      unsigned int eidx2 = idx2*nelt_;
//      d_m2_idx = (double)eidx2;
//      mat = Mat_->begin() + (int)(0.5*d_m2_idx*(TwoN-d_m2_idx-1.0)+d_m2_idx);
      double Vj = vect_[eidx2 + jidx];
      for (unsigned int idx1 = idx2; idx1 < mass_.size(); idx1++) {
        double Mass = sqrt( mass2 * mass_[idx1] );
        if (idx2 == idx1) {
          // Self
          for (unsigned int ej = jidx; ej < nelt_; ej++) {
            *mat = (*mat - (Vj * vect_[eidx2 + ej])) * Mass;
            ++mat;
          }
        } else {
          unsigned int eidx1 = idx1*nelt_;
          for (unsigned int iidx = 0; iidx < nelt_; iidx++) {
            *mat = (*mat - (Vj * vect_[eidx1 + iidx])) * Mass;
            ++mat;
          }
        }
      } // END inner loop over idx1
    } // END loop over elements of vect_[idx2]
  } // END outer loop over idx2
  return 0;
}

#include <algorithm>
#include "InternalCoords.h"
//#incl ude "../Vec3.h"

using namespace Cpptraj::Structure;

const int InternalCoords::NO_ATOM = -1;

/** CONSTRUCTOR */
InternalCoords::InternalCoords() {
  std::fill(idx_, idx_+3, NO_ATOM);
  std::fill(val_, val_+3, 0);
//  std::fill(xyz_, xyz_+3, 0);
}

/** CONSTRUCTOR - Take indices and values */
InternalCoords::InternalCoords(int atJ, int atK, int atL, double dist, double theta, double phi) {
  if (atJ < 0)
    idx_[0] = NO_ATOM;
  else
    idx_[0] = atJ;
  if (atK < 0)
    idx_[1] = NO_ATOM;
  else
    idx_[1] = atK;
  if (atL < 0)
    idx_[2] = NO_ATOM;
  else
    idx_[2] = atL;

  val_[0] = dist;
  val_[1] = theta;
  val_[2] = phi;
}

/** CONSTRUCTOR - pointer to XYZ coords, for first seed atom */
/*InternalCoords::InternalCoords(const double* xyz) {
  std::fill(idx_, idx_+3, NO_ATOM);
  std::fill(val_, val_+3, 0);
//  std::copy( xyz, xyz+3, xyz_ );
//  isSet_ = true;
}*/

/** COPY CONSTRUCTOR */
InternalCoords::InternalCoords(InternalCoords const& rhs) {
  std::copy( rhs.idx_, rhs.idx_+3, idx_ );
  std::copy( rhs.val_, rhs.val_+3, val_ );
//  std::copy( rhs.xyz_, rhs.xyz_+3, xyz_ );
}

/** ASSIGNMENT */
InternalCoords& InternalCoords::operator=(InternalCoords const& rhs) {
  if (&rhs == this) return *this;
  std::copy( rhs.idx_, rhs.idx_+3, idx_ );
  std::copy( rhs.val_, rhs.val_+3, val_ );
  //std::copy( rhs.xyz_, rhs.xyz_+3, xyz_ );
  //isSet_ = rhs.isSet_;
  return *this;
}

/** Zero out the XYZ coordinates. */
/*void InternalCoords::ZeroXYZ() {
  std::fill(idx_, idx_+3, NO_ATOM);
  std::fill(val_, val_+3, 0);
  std::fill(xyz_, xyz_+3, 0);
  isSet_ = true;
}*/

/** Set xyz coords */
/*void InternalCoords::SetXYZ(Vec3 const& xyz) {
  xyz_[0] = xyz[0];
  xyz_[1] = xyz[1];
  xyz_[2] = xyz[2];
  isSet_ = true;
}*/

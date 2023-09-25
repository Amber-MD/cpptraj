#include <algorithm>
#include "InternalCoords.h"

using namespace Cpptraj::Structure;

const int InternalCoords::NO_ATOM = -1;

/** CONSTRUCTOR */
InternalCoords::InternalCoords() : isSet_(false) {
  std::fill(idx_, idx_+3, NO_ATOM);
  std::fill(val_, val_+3, 0);
  std::fill(xyz_, xyz_+3, 0);
}

/** CONSTRUCTOR - pointer to XYZ coords, for first seed atom */
InternalCoords::InternalCoords(const double* xyz) {
  std::fill(idx_, idx_+3, NO_ATOM);
  std::fill(val_, val_+3, 0);
  std::copy( xyz, xyz+3, xyz_ );
  isSet_ = true;
}

/** COPY CONSTRUCTOR */
InternalCoords::InternalCoords(InternalCoords const& rhs) : isSet_(rhs.isSet_) {
  std::copy( rhs.idx_, rhs.idx_+3, idx_ );
  std::copy( rhs.val_, rhs.val_+3, val_ );
  std::copy( rhs.xyz_, rhs.xyz_+3, xyz_ );
}

/** ASSIGNMENT */
InternalCoords& InternalCoords::operator=(InternalCoords const& rhs) {
  if (&rhs == this) return *this;
  std::copy( rhs.idx_, rhs.idx_+3, idx_ );
  std::copy( rhs.val_, rhs.val_+3, val_ );
  std::copy( rhs.xyz_, rhs.xyz_+3, xyz_ );
  isSet_ = rhs.isSet_;
  return *this;
}


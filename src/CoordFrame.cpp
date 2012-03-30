#include "CoordFrame.h"

// CONSTRUCTOR
CoordFrame::CoordFrame() {
}

// CONSTRUCTOR - double* and natoms
CoordFrame::CoordFrame(int natom, const double *Xin) {
  coords_.reserve( natom );
  // NOTE: No NULL CHECK!
  double *Xptr = (double*)Xin;
  for (int atom = 0; atom < natom; ++atom) {
    coords_.push_back( Xptr );
    Xptr += 3;
  }
}

// COPY CONSTRUCTOR
CoordFrame::CoordFrame(const CoordFrame &rhs) :
  coords_(rhs.coords_)
{}

// Assignment Operator
CoordFrame &CoordFrame::operator=(const CoordFrame &rhs) {
  if (&rhs == this) return *this;
  coords_ = rhs.coords_;
  return *this;
}

double CoordFrame::DIST2_NoImage(int a1, int a2) {
  // NOTE: no bounds check!
  Vec3 delta = coords_[a1];
  delta -= coords_[a2];
  return delta.LengthSquared();
}


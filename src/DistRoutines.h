#ifndef INC_DISTROUTINES_H
#define INC_DISTROUTINES_H
class Matrix_3x3;
#include "Vec3.h"
#include "Box.h"
#include "ImageOption.h"
#include <cmath> // sqrt()
/*! \file DistRoutines.h
    \brief A collection of routines used to calculate distance.
 */
/// \return Vector representing minimum imaged distance between two points.
Vec3 MinImagedVec(Vec3 const&, Vec3 const&, Matrix_3x3 const&, Matrix_3x3 const&);

/// \return Coordinates of minimum image of a2 to a1
Vec3 MinImagedCoords(double&, ImageOption::Type, Vec3 const&, Vec3 const&, Box const&);

/// \return Distance squared between points in fractional space, general unit cell imaging.
double DIST2_ImageNonOrthoRecip(Vec3 const&, Vec3 const&, double, int*, Matrix_3x3 const&);
/// \return Distance squared between points in Cartesian space, general unit cell imaging.
double DIST2_ImageNonOrtho(Vec3 const&, Vec3 const&, Matrix_3x3 const&, Matrix_3x3 const&);
/// \return Distance squared between points in Cartesian space, general unit cell imaging.
double DIST2_ImageNonOrtho(const double*, const double*, Matrix_3x3 const&, Matrix_3x3 const&);

// Frame::DIST2_ImageOrtho()
/** Return the minimum orthorhombic imaged distance^2 between coordinates a1 
  * and a2.
  */
inline double DIST2_ImageOrtho(const double* a1, const double* a2, Box const& box)
{
  double boxx = box.Param(Box::X);
  double boxy = box.Param(Box::Y);
  double boxz = box.Param(Box::Z);
  // If box lengths are zero no imaging possible
  if (boxx==0.0 || boxy==0.0 || boxz==0.0) return -1.0;
  double x = a1[0] - a2[0];
  double y = a1[1] - a2[1];
  double z = a1[2] - a2[2];
  // Get rid of sign info
  if (x<0) x=-x;
  if (y<0) y=-y;
  if (z<0) z=-z;
  // Get rid of multiples of box lengths 
  while (x > boxx) x = x - boxx;
  while (y > boxy) y = y - boxy;
  while (z > boxz) z = z - boxz;
  // Find shortest distance in periodic reference
  double D = boxx - x;
  if (D < x) x = D;
  D = boxy - y;
  if (D < y) y = D;  
  D = boxz - z;
  if (D < z) z = D;

  return (x*x + y*y + z*z);
}

/** Return the minimum orthorhombic imaged distance^2 between coordinates a1 
  * and a2.
  */
inline double DIST2_ImageOrtho(Vec3 const& a1, Vec3 const& a2, Box const& box) {
  return DIST2_ImageOrtho(a1.Dptr(), a2.Dptr(), box);
}

// Frame::DIST2_NoImage()
/** Return distance^2 between coordinates in a1 and a2.
  */
inline double DIST2_NoImage(const double* a1, const double* a2) {
  double x = a1[0] - a2[0];
  double y = a1[1] - a2[1];
  double z = a1[2] - a2[2];
  //double D = x*x + y*y + z*z;
  //fprintf(stdout,"Mask1=%8.3f %8.3f %8.3f Mask2=%8.3f %8.3f %8.3f D=%8.3f\n",
  //        a1[0],a1[1],a1[2],a2[0],a2[1],a2[2],D);
  return (x*x + y*y + z*z);
}

/// \return Distance squared, no imaging.
inline double DIST2_NoImage( Vec3 const& a1, Vec3 const& a2 ) {
  Vec3 vec = a1 - a2;
  return vec.Magnitude2();
}

/// \return Distance, no imaging.
inline double DIST_NoImage( Vec3 const& a1, Vec3 const& a2 ) {
  Vec3 vec = a1 - a2;
  return sqrt( vec.Magnitude2() );
}

/** \return Distance squared using either minimum-image convention or no imaging. */
inline double DIST2(ImageOption::Type itype, const double* a1, const double* a2, Box const& box) {
/*  if (imagingEnabled) {
    if (is_ortho)
      return DIST2_ImageOrtho(a1, a2, box);
    else
      return DIST2_ImageNonOrtho(a1, a2, box.UnitCell(), box.FracCell());
  } else
    return DIST2_NoImage(a1, a2);*/
  if (itype == ImageOption::NO_IMAGE)
    return DIST2_NoImage(a1, a2);
  else if (itype == ImageOption::ORTHO)
    return DIST2_ImageOrtho(a1, a2, box);
  else // NONORTHO
    return DIST2_ImageNonOrtho(a1, a2, box.UnitCell(), box.FracCell());
}

/** \return Distance squared using either minimum-image convention or no imaging. */
inline double DIST2(ImageOption::Type itype, Vec3 const& a1, Vec3 const& a2, Box const& box) {
  return DIST2(itype, a1.Dptr(), a2.Dptr(), box);
}

/** \return Distance using either minimum-image convention or no imaging. */
inline double DIST(ImageOption::Type itype, Vec3 const& a1, Vec3 const& a2, Box const& box) {
  return sqrt(DIST2(itype, a1.Dptr(), a2.Dptr(), box));
}
#endif

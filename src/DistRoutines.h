#ifndef INC_DISTROUTINES_H
#define INC_DISTROUTINES_H
#include <cmath>
#include "Box.h"
#include "ImageOption.h"
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

/// \return Distance squared, X-aligned and orthogonal imaging.
inline double DIST2_ImageOrtho(const double* a1, const double* a2,
                               Box const& box)
{
  const double boxX = box.Param(Box::X);
  const double boxY = box.Param(Box::Y);
  const double boxZ = box.Param(Box::Z);

  if (boxX == 0.0 || boxY == 0.0 || boxZ == 0.0)
    return -1.0;

  double x = std::fabs(a1[0] - a2[0]);
  double y = std::fabs(a1[1] - a2[1]);
  double z = std::fabs(a1[2] - a2[2]);

  while (x > boxX) x -= boxX;
  while (y > boxY) y -= boxY;
  while (z > boxZ) z -= boxZ;

  double D = boxX - x;
  if (D < x) x = D;
  D = boxY - y;
  if (D < y) y = D;
  D = boxZ - z;
  if (D < z) z = D;

  return x*x + y*y + z*z;
}
/// \return Distance squared, X-aligned and orthogonal imaging.
double DIST2_ImageOrtho(Vec3 const&, Vec3 const&, Box const&);

/// \return Distance squared, no imaging.
inline double DIST2_NoImage(const double* a1, const double* a2)
{
  const double x = a1[0] - a2[0];
  const double y = a1[1] - a2[1];
  const double z = a1[2] - a2[2];
  return x*x + y*y + z*z;
}
/// \return Distance squared, no imaging
double DIST2_NoImage( Vec3 const&, Vec3 const& );
/// \return Distance, no imaging.
double DIST_NoImage( Vec3 const&, Vec3 const& );

/// \return Distance squared using minimum-image convention or no imaging.
inline double DIST2(ImageOption::Type itype, const double* a1,
                    const double* a2, Box const& box)
{
  if (itype == ImageOption::NO_IMAGE)
    return DIST2_NoImage(a1, a2);
  else if (itype == ImageOption::ORTHO)
    return DIST2_ImageOrtho(a1, a2, box);
  else
    return DIST2_ImageNonOrtho(a1, a2, box.UnitCell(), box.FracCell());
}
// \return Distance squared using minimum-image convention or no imaging.
double DIST2(ImageOption::Type, Vec3 const&, Vec3 const&, Box const&);
// \return Distance using minimum-image convention or no imaging.
double DIST(ImageOption::Type, Vec3 const&, Vec3 const&, Box const&);
#endif

#ifndef INC_DISTROUTINES_H
#define INC_DISTROUTINES_H
class Vec3;
class Matrix_3x3;
class Box;
/*! \file DistRoutines.h
    \brief A collection of routines used to calculate distance.
 */
// TODO: Deprecate ImagingType, just use Box::Is_X_Aligned_Ortho() 
/// Potential imaging types 
enum ImagingType { NOIMAGE=0, ORTHO, NONORTHO };
/// \return Vector representing minimum imaged distance between two points.
Vec3 MinImagedVec(Vec3 const&, Vec3 const&, Matrix_3x3 const&, Matrix_3x3 const&);
/// \return Distance squared between points in Cartesian space, general unit cell imaging.
double DIST2_ImageNonOrtho(Vec3 const&, Vec3 const&, Matrix_3x3 const&, Matrix_3x3 const&);
/// \return Distance squared between points in fractional space, general unit cell imaging.
double DIST2_ImageNonOrthoRecip(Vec3 const&, Vec3 const&, double, int*, Matrix_3x3 const&);
/// \return Distance squared, X-aligned and orthogonal imaging.
double DIST2_ImageOrtho(Vec3 const&, Vec3 const&, Box const&);
/// \return Distance squared, no imaging.
double DIST2_NoImage(const double*, const double*);
/// \return Distance squared, no imaging
double DIST2_NoImage( Vec3 const&, Vec3 const& );
/// \return Distance, no imaging.
double DIST_NoImage( Vec3 const&, Vec3 const& );
/// \return Distance squared, optionally imaged.
double DIST2(const double*, const double*, ImagingType, Box const&, // TODO deprecate this version 
             Matrix_3x3 const&, Matrix_3x3 const&);
/// \return Distance squared using minimum-image convention or no imaging.
double DIST2(bool, Vec3 const&, Vec3 const&, Box const&);
/// \return Distance using minimum-image convention or no imaging.
double DIST(bool, Vec3 const&, Vec3 const&, Box const&);
#endif

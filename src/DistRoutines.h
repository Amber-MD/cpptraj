#ifndef INC_DISTROUTINES_H
#define INC_DISTROUTINES_H
#include "Vec3.h"
/*! \file DistRoutines.h
    \brief A collection of routines used to calculate distance.
 */
/// Potential imaging types 
enum ImagingType { NOIMAGE=0, ORTHO, NONORTHO };
double DIST2_ImageNonOrtho(const double *, const double *, const double *, const double *);
double DIST2_ImageNonOrthoRecip(const double *, const double *, double minIn, 
                                int *, const double *);
double DIST2_ImageOrtho(const double *, const double *, Vec3 const&);
double DIST2_NoImage(const double *, const double *);
double DIST2(const double*, const double*, ImagingType, Vec3 const&, 
             const double*, const double*);
#endif

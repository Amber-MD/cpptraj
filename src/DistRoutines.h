#ifndef INC_DISTROUTINES_H
#define INC_DISTROUTINES_H
/*! \file DistRoutines.h
    \brief A collection of routines used to calculate distance.
 */
double MinImageNonOrtho2(double *Coord1, double *Coord2, double *box, int origin, int *ixyz,
                         double *ucell, double *recip);
double DIST2_ImageNonOrtho(double *a1, double *a2, double *ucell, double *recip);
double DIST2_ImageNonOrthoRecip(double *f, double *f2, double minIn, int *ixyz, double *ucell);
double DIST2_ImageOrtho(double *a1, double *a2, double *box);
double DIST2_NoImage(const double *a1, const double *a2);

#endif

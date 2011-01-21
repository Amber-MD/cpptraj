#ifndef INC_DISTROUTINES_H
#define INC_DISTROUTINES_H
// DistRoutines
double MinImageNonOrtho2(double *Coord1, double *Coord2, double *box, bool origin, int *ixyz,
                         double *ucell, double *recip);
double DIST2_ImageNonOrtho(double *a1, double *a2, double *ucell, double *recip);
double DIST2_ImageNonOrtho(double *f, double *f2, double minIn, int *ixyz, double *ucell);
double DIST2_ImageOrtho(double *a1, double *a2, double *box);
double DIST2_NoImage(double *a1, double *a2);
#endif

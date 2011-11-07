#ifndef INC_DISTROUTINES_H
#define INC_DISTROUTINES_H
// DistRoutines
double MinImageNonOrtho2(double *Coord1, double *Coord2, double *box, int origin, int *ixyz,
                         double *ucell, double *recip);
#ifdef __cplusplus
extern "C"
#endif
double DIST2_ImageNonOrtho(double *a1, double *a2, double *ucell, double *recip);

double DIST2_ImageNonOrthoRecip(double *f, double *f2, double minIn, int *ixyz, double *ucell);

#ifdef __cplusplus
extern "C"
#endif
double DIST2_ImageOrtho(double *a1, double *a2, double *box);

#ifdef __cplusplus
extern "C"
#endif
double DIST2_NoImage(double *a1, double *a2);

#ifdef __cplusplus
extern "C" 
#endif
double boxToRecip(double *, double *, double *);

#endif

#ifndef INC_TORSIONROUTINES_H
#define INC_TORSIONROUTINES_H
/// TorsionRoutines
/// A collection of routines used to calculate torsions, adapted from Ptraj.
#ifdef __cplusplus
extern "C"
#endif
double Torsion(double *, double *, double *, double *);

double Pucker_AS(double *, double *, double *, double *, double *, double *);
double Pucker_CP(double *, double *, double *, double *, double *, double *);
double CalcAngle(double*, double*, double*);
#endif

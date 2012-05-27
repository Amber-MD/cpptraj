#ifndef INC_TORSIONROUTINES_H
#define INC_TORSIONROUTINES_H
/*! \file TorsionRoutines.h
    \brief A collection of routines used to calculate angles and torsions
 */

double Torsion(double *, double *, double *, double *);

double Pucker_AS(double *, double *, double *, double *, double *, double *);
double Pucker_CP(double *, double *, double *, double *, double *, double *);
double CalcAngle(double*, double*, double*);
#endif

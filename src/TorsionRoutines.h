#ifndef INC_TORSIONROUTINES_H
#define INC_TORSIONROUTINES_H
/*! \file TorsionRoutines.h
    \brief A collection of routines used to calculate angles and torsions
 */
class Vec3;
double Torsion(const double *, const double *, const double *, const double *);

void Torsion_and_part_deriv(const double*, const double*, const double*, const double*,
                            Vec3&, Vec3&, Vec3&, Vec3&, double&, double&);

double Pucker_AS(const double*, const double*, const double*, const double*, 
                 const double*, double&);
double Pucker_CP(const double*, const double*, const double*, const double*, 
                 const double*, const double*, int, double&, double&);
double CalcAngle(const double*, const double*, const double*);
#endif

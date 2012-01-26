/*! \file TorsionRoutines.cpp
    \brief Routines used to calculate torsions and angles.
 */
#include <cmath>
#include "vectormath.h"
#include "TorsionRoutines.h"
#include "Constants.h" // PI, TWOPI

// Torsion()
/** Given 4 sets of XYZ coords, calculate the torsion (in radians) between the 
  * planes formed by a1-a2-a3 and a2-a3-a4.
  */
extern "C" double Torsion(double *a1, double *a2, double *a3, double *a4) {
  double Lx, Ly, Lz, Lnorm;
  double Rx, Ry, Rz, Rnorm;
  double Sx, Sy, Sz;
  double angle;

  CROSS_PRODUCT(     Lx,      Ly,      Lz,
                (a2[0]-a1[0]), (a2[1]-a1[1]), (a2[2]-a1[2]),
                (a3[0]-a2[0]), (a3[1]-a2[1]), (a3[2]-a2[2]));

  CROSS_PRODUCT(     Rx,      Ry,      Rz,
                (a4[0]-a3[0]), (a4[1]-a3[1]), (a4[2]-a3[2]),
                (a2[0]-a3[0]), (a2[1]-a3[1]), (a2[2]-a3[2]));

  Lnorm = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
  Rnorm = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);

  CROSS_PRODUCT(Sx, Sy, Sz,
                Lx, Ly, Lz,
                Rx, Ry, Rz);

  angle = (Lx*Rx + Ly*Ry + Lz*Rz) / (Lnorm * Rnorm);

  if ( angle > 1.0 ) angle = 1.0;
  if ( angle < -1.0 ) angle = -1.0;

  angle = acos( angle );

  if ( (Sx * (a3[0]-a2[0]) + Sy * (a3[1]-a2[1]) + Sz * (a3[2]-a2[2])) < 0 )
    angle = -angle;

  return angle;
}

// Pucker_AS()
/** Return the pucker (in radians) of coords stored in a1-a5 based on 
  * Altona & Sundaralingam method.
  */
double Pucker_AS(double *a1, double *a2, double *a3, double *a4, double *a5, double *amp) {
  double pucker;
  double v1, v2, v3, v4, v5, a, b;
  double pi_over_5;

  pucker = 0.0;
  *amp = 0.0;
  pi_over_5 = PI / 5.0;

  v4 = Torsion(a4,a5,a1,a2);
  v5 = Torsion(a5,a1,a2,a3);
  v1 = Torsion(a1,a2,a3,a4);
  v2 = Torsion(a2,a3,a4,a5);
  v3 = Torsion(a3,a4,a5,a1);

  a = (v1*cos(0.0) +
       v2*cos( 4.0*pi_over_5) +
       v3*cos( 8.0*pi_over_5) +
       v4*cos(12.0*pi_over_5) +
       v5*cos(16.0*pi_over_5))*0.4;

  b = (v1*sin(0.0) +
       v2*sin( 4.0*pi_over_5) +
       v3*sin( 8.0*pi_over_5) +
       v4*sin(12.0*pi_over_5) +
       v5*sin(16.0*pi_over_5))*-0.4;

  *amp = sqrt(a*a + b*b);

  if (*amp != 0.0)
    pucker = atan2(b,a);
  if (pucker < 0) pucker += TWOPI;

  return pucker;
}

// Pucker_CP()
/** Return the pucker (in radians) of coords in a1-a5 based on method of
  * Cremer & Pople.
  */
double Pucker_CP(double *a1, double *a2, double *a3, double *a4, double *a5, double *amplitude) {
  double pucker, norm;
  double x1, y1, z1;
  double x2, y2, z2;
  double x3, y3, z3;
  double x4, y4, z4;
  double x5, y5, z5;
  double nx, ny, nz;
  double rcx, rcy, rcz;
  double r1x, r1y, r1z;
  double r2x, r2y, r2z;
  double sum1, sum2;
  double pi_over_5;
  double one_over_five;

  pucker = 0.0;
  *amplitude = 0.0;
  one_over_five = 1 / 5.0;
  pi_over_5 = PI * one_over_five;

  x2 = a1[0]; y2 = a1[1]; z2 = a1[2];
  x3 = a2[0]; y3 = a2[1]; z3 = a2[2];
  x4 = a3[0]; y4 = a3[1]; z4 = a3[2];
  x5 = a4[0]; y5 = a4[1]; z5 = a4[2];
  x1 = a5[0]; y1 = a5[1]; z1 = a5[2];

  // Calculate geometric center
  rcx = (x1 + x2 + x3 + x4 + x5)*one_over_five;
  rcy = (y1 + y2 + y3 + y4 + y5)*one_over_five;
  rcz = (z1 + z2 + z3 + z4 + z5)*one_over_five;

  x1 -= rcx; y1 -= rcy; z1 -=rcz;
  x2 -= rcx; y2 -= rcy; z2 -=rcz;
  x3 -= rcx; y3 -= rcy; z3 -=rcz;
  x4 -= rcx; y4 -= rcy; z4 -=rcz;
  x5 -= rcx; y5 -= rcy; z5 -=rcz;

  // Calculate normal vectors
  r1x = x1 * sin(0.0) +
        x2 * sin(2.0*pi_over_5) +
        x3 * sin(4.0*pi_over_5) +
        x4 * sin(6.0*pi_over_5) +
        x5 * sin(8.0*pi_over_5);
  r1y = y1 * sin(0.0) +
        y2 * sin(2.0*pi_over_5) +
        y3 * sin(4.0*pi_over_5) +
        y4 * sin(6.0*pi_over_5) +
        y5 * sin(8.0*pi_over_5);
  r1z = z1 * sin(0.0) +
        z2 * sin(2.0*pi_over_5) +
        z3 * sin(4.0*pi_over_5) +
        z4 * sin(6.0*pi_over_5) +
        z5 * sin(8.0*pi_over_5);

  r2x = x1 * cos(0.0) +
        x2 * cos(2.0*pi_over_5) +
        x3 * cos(4.0*pi_over_5) +
        x4 * cos(6.0*pi_over_5) +
        x5 * cos(8.0*pi_over_5);
  r2y = y1 * cos(0.0) +
        y2 * cos(2.0*pi_over_5) +
        y3 * cos(4.0*pi_over_5) +
        y4 * cos(6.0*pi_over_5) +
        y5 * cos(8.0*pi_over_5);
  r2z = z1 * cos(0.0) +
        z2 * cos(2.0*pi_over_5) +
        z3 * cos(4.0*pi_over_5) +
        z4 * cos(6.0*pi_over_5) +
        z5 * cos(8.0*pi_over_5);

  // Calculate vector normal to plane
  CROSS_PRODUCT( nx,  ny,  nz,
                r1x, r1y, r1z,
                r2x, r2y, r2z );

  norm = sqrt(nx*nx + ny*ny + nz*nz);
  nx /= norm;
  ny /= norm;
  nz /= norm;

  // Rotate around Z axis
  z1 = x1*nx + y1*ny + z1*nz;
  z2 = x2*nx + y2*ny + z2*nz;
  z3 = x3*nx + y3*ny + z3*nz;
  z4 = x4*nx + y4*ny + z4*nz;
  z5 = x5*nx + y5*ny + z5*nz;

  sum1 = z1 * cos(0.0) +
         z2 * cos( 4.0*pi_over_5) +
         z3 * cos( 8.0*pi_over_5) +
         z4 * cos(12.0*pi_over_5) +
         z5 * cos(16.0*pi_over_5);
  sum2 = -(z1 * sin(0.0) +
           z2 * sin( 4.0*pi_over_5) +
           z3 * sin( 8.0*pi_over_5) +
           z4 * sin(12.0*pi_over_5) +
           z5 * sin(16.0*pi_over_5));

  norm = sqrt(sum1*sum1 + sum2*sum2);
  *amplitude = norm * sqrt(2.0*one_over_five);
  pucker = asin( sum2 / norm );
  if (sum1 < 0.0)
    pucker = PI - pucker;
  else if (pucker < 0.0)
    pucker += TWOPI;

  return pucker;
}

// CalcAngle()
double CalcAngle(double *V1, double *V2, double *V3)
{
  double angle;
  double xij = V1[0] - V2[0];
  double yij = V1[1] - V2[1];
  double zij = V1[2] - V2[2];
  
  double xkj = V3[0] - V2[0];
  double ykj = V3[1] - V2[1];
  double zkj = V3[2] - V2[2];
  
  double rij = xij*xij + yij*yij + zij*zij;
  double rkj = xkj*xkj + ykj*ykj + zkj*zkj;

  if (rij > SMALL && rkj > SMALL) {
    angle = (xij*xkj + yij*ykj + zij*zkj) / sqrt(rij*rkj);
    if (angle > 1.0)
      angle = 1.0;
    else if (angle < -1.0)
      angle = -1.0;
    angle = acos(angle);
  } else
    angle = 0.0;

  return angle;
}


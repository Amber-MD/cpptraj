/*! \file TorsionRoutines.cpp
    \brief Routines used to calculate torsions and angles.
 */
#include <cmath>
#include "TorsionRoutines.h"
#include "Constants.h" // PI, TWOPI
#include "Vec3.h"
#ifdef NEW_PUCKER_CODE
#  include "CpptrajStdio.h" // DEBUG
#endif

// Torsion()
/** Given 4 sets of XYZ coords, calculate the torsion (in radians) between the 
  * planes formed by a1-a2-a3 and a2-a3-a4.
  */
double Torsion(const double *a1, const double *a2, const double *a3, const double *a4) 
{
  double Lx = ((a2[1]-a1[1])*(a3[2]-a2[2])) - ((a2[2]-a1[2])*(a3[1]-a2[1])); 
  double Ly = ((a2[2]-a1[2])*(a3[0]-a2[0])) - ((a2[0]-a1[0])*(a3[2]-a2[2])); 
  double Lz = ((a2[0]-a1[0])*(a3[1]-a2[1])) - ((a2[1]-a1[1])*(a3[0]-a2[0]));

  double Rx = ((a4[1]-a3[1])*(a2[2]-a3[2])) - ((a4[2]-a3[2])*(a2[1]-a3[1])); 
  double Ry = ((a4[2]-a3[2])*(a2[0]-a3[0])) - ((a4[0]-a3[0])*(a2[2]-a3[2])); 
  double Rz = ((a4[0]-a3[0])*(a2[1]-a3[1])) - ((a4[1]-a3[1])*(a2[0]-a3[0]));

  double Lnorm = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
  double Rnorm = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);

  double Sx = (Ly*Rz) - (Lz*Ry); 
  double Sy = (Lz*Rx) - (Lx*Rz); 
  double Sz = (Lx*Ry) - (Ly*Rx);

  double angle = (Lx*Rx + Ly*Ry + Lz*Rz) / (Lnorm * Rnorm);

  if ( angle > 1.0 ) angle = 1.0;
  if ( angle < -1.0 ) angle = -1.0;

  angle = acos( angle );

  if ( (Sx * (a3[0]-a2[0]) + Sy * (a3[1]-a2[1]) + Sz * (a3[2]-a2[2])) < 0 )
    angle = -angle;

  return angle;
}

/// Constant used in CP pucker calc
static const double one_over_five = 1.0 / 5.0;
/// Constant used in AS and CP pucker calcs
static const double pi_over_5 = PI * one_over_five;

// Pucker_AS()
/** Return the pucker (in radians) of coords stored in a1-a5 based on 
  * Altona & Sundaralingam method.
  */
double Pucker_AS(const double* a1, const double* a2, const double* a3, 
                 const double* a4, const double* a5, double& amp) 
{
  double pucker;
  double v1, v2, v3, v4, v5, a, b;

  pucker = 0.0;
  amp = 0.0;

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

  amp = sqrt(a*a + b*b);

  if (amp != 0.0)
    pucker = atan2(b,a);
  if (pucker < 0) pucker += TWOPI;

  return pucker;
}

// Pucker_CP()
/** Return the pucker (in radians) of coords in a1-a5 based on method of
  * Cremer & Pople.
  */
double Pucker_CP(const double* a1, const double* a2, const double* a3, 
                 const double* a4, const double* a5, const double* a6,
                 int N, double& amplitude) 
{
# ifndef NEW_PUCKER_CODE
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

  pucker = 0.0;
  amplitude = 0.0;
  x2 = a1[0]; y2 = a1[1]; z2 = a1[2];
  x3 = a2[0]; y3 = a2[1]; z3 = a2[2];
  x4 = a3[0]; y4 = a3[1]; z4 = a3[2];
  x5 = a4[0]; y5 = a4[1]; z5 = a4[2];
  x1 = a5[0]; y1 = a5[1]; z1 = a5[2];

  // Calculate geometric center
  rcx = (x1 + x2 + x3 + x4 + x5)*one_over_five;
  rcy = (y1 + y2 + y3 + y4 + y5)*one_over_five;
  rcz = (z1 + z2 + z3 + z4 + z5)*one_over_five;
  // Translate to center
  x1 -= rcx; y1 -= rcy; z1 -=rcz;
  x2 -= rcx; y2 -= rcy; z2 -=rcz;
  x3 -= rcx; y3 -= rcy; z3 -=rcz;
  x4 -= rcx; y4 -= rcy; z4 -=rcz;
  x5 -= rcx; y5 -= rcy; z5 -=rcz;
  // Calculate position vectors
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
  nx = (r1y*r2z) - (r1z*r2y); 
  ny = (r1z*r2x) - (r1x*r2z); 
  nz = (r1x*r2y) - (r1y*r2x);
  // Normalize
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
  amplitude = norm * sqrt(2.0*one_over_five);
  pucker = asin( sum2 / norm );
# else
// -------------------------------------
  double dN = (double)N;
  double twopi_over_N = TWOPI / dN;
  Vec3 XYZ[6];
  XYZ[1].Assign( a1 ); 
  XYZ[2].Assign( a2 ); 
  XYZ[3].Assign( a3 ); 
  XYZ[4].Assign( a4 );
  if (N == 5) {
    XYZ[0].Assign( a5 ); // Ring apex
  } else if (N == 6) {
    XYZ[0].Assign( a6 ); // Ring apex
    XYZ[5].Assign( a5 );
  } else
    return -1.0; // Internal error
  // Calculate geometric center
  Vec3 Center(0.0, 0.0, 0.0);
  for (int i = 0; i < N; i++)
    Center += XYZ[i];
  Center /= dN;
  //Center.Print("RCXYZ");
  // Translate to center
  for (int i = 0; i < N; i++)
    XYZ[i] -= Center;
  //XYZ[0].Print("Trans[0]");
  //XYZ[1].Print("Trans[1]");
  //XYZ[2].Print("Trans[2]");
  //XYZ[3].Print("Trans[3]");
  //XYZ[4].Print("Trans[4]");
  // Calculate position vectors
  Vec3 R1(0.0, 0.0, 0.0);
  Vec3 R2(0.0, 0.0, 0.0);
  for (int i = 0; i < N; i++) {
    double factor = twopi_over_N * (double)i;
    double sin_val = sin( factor );
    double cos_val = cos( factor );
    R1[0] += XYZ[i][0] * sin_val;
    R2[0] += XYZ[i][0] * cos_val;
    R1[1] += XYZ[i][1] * sin_val;
    R2[1] += XYZ[i][1] * cos_val;
    R1[2] += XYZ[i][2] * sin_val; 
    R2[2] += XYZ[i][2] * cos_val;
  }
  //R1.Print("R1");
  //R2.Print("R2");
  // Calculate vector normal to plane
  Vec3 NXYZ = R1.Cross( R2 );
  // Normalize
  NXYZ.Normalize();
  //NXYZ.Print("normN");
  // Rotate around Z axis
  double sum1 = 0.0;
  double sum2 = 0.0;
  double sum3 = 0.0;
  double Zn[6];
  double sumZ = 0.0; // DEBUG
  //mprintf("Z=");
  for (int i = 0; i < N; i++) {
    double factor = 2.0 * twopi_over_N * (double)i;
    double cos_val = cos( factor );
    double sin_val = sin( factor ); 
    Zn[i] = XYZ[i] * NXYZ;
    sumZ += Zn[i] * Zn[i]; // DEBUG
    mprintf("DEBUG:\t Z[%i]= %f\n", i+1, Zn[i]);
    sum1 += Zn[i] * cos_val;
    sum2 -= Zn[i] * sin_val;
  }
  mprintf("DEBUG: sqrt(sumZ^2)= %f\n", sqrt(sumZ));
  // For even # coords (only 6 currently) calc extra pucker coord
  if (N == 6) {
    double mult = 1.0;
    for (int i = 0; i < N; i++) {
      sum3 += mult * Zn[i]; // mult ~ pow( -1.0, i )
      mult = -mult;
    }
    sum3 /= sqrt( dN );
    mprintf("DEBUG: sum3= %f\n", sum3);
  }
  //mprintf("\n");
  mprintf("DEBUG: sum1= %f   sum2= %f\n", sum1, sum2);
  double norm = sqrt(sum1*sum1 + sum2*sum2);
  amplitude = norm * sqrt( 2.0 / dN );
  double pucker = asin( sum2 / norm );
  // DEBUG - Recover Zn values
  for (int i = 0; i < N; i++) {
    double Zi = sqrt(2.0/dN) * ((sqrt(2.0/dN)*sum1)/cos(pucker)) * 
                cos(pucker + ((FOURPI * (double)i)/dN));
    if (N == 6) Zi += (1.0/sqrt(dN)) * sum3 * pow(-1.0, i);
    mprintf("DEBUG:\tCalcd. Zn[%i]= %f\n", i+1, Zi);
  }
# endif    
  if (sum1 < 0.0)
    pucker = PI - pucker;
  else if (pucker < 0.0)
    pucker += TWOPI;

  return pucker;
}

// CalcAngle()
double CalcAngle(const double* V1, const double* V2, const double* V3)
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

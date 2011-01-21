#include <cmath>
#include "DistRoutines.h"
#include "vectormath.h" // DEGRAD

/*
 * Frame::ClosestImage()
 * Given two coordinates A and B, determine the unit XYZ vector that points 
 * towards the closest image of B to A.
 * It is assumed the coordinates are already relative to the box center.
 */
/*
void Frame::ClosestImage(double *A, double *B, int *ixyz) {
  double halfBox[3];//, delta;
  int vectorA[3], vectorB[3], i;

  mprintf("DEBUG: CoordA     = %lf %lf %lf\n",A[0],A[1],A[2]);
  mprintf("DEBUG: CoordB     = %lf %lf %lf\n",B[0],B[1],B[2]);

  halfBox[0] = box[0] / 6.0; 
  halfBox[1] = box[1] / 6.0; 
  halfBox[2] = box[2] / 6.0;
  mprintf("DEBUG: Half box = %lf %lf %lf\n",halfBox[0],halfBox[1],halfBox[2]);

  // Vector A
  vectorA[0] = 0; vectorA[1] = 0; vectorA[2] = 0;
  for (i=0; i<3; i++) {
    if (A[i] < -halfBox[i]) vectorA[i] = -1;
    if (A[i] >  halfBox[i]) vectorA[i] =  1;
//    delta = A[i] - boxCenter[i];
//  if (delta > 0.0) vectorA[i] = 1;
//  if (delta < 0.0) vectorA[i] = -1;
//
  }
  mprintf("DEBUG:  VectorA = %2i %2i %2i\n",vectorA[0],vectorA[1],vectorA[2]);

  // NOT Vector B
  vectorB[0] = 0; vectorB[1] = 0; vectorB[2] = 0;
  for (i=0; i<3; i++) {
    if (B[i] < -halfBox[i]) vectorB[i] =  1; // NOT
    if (B[i] >  halfBox[i]) vectorB[i] = -1; // NOT
//    delta = B[i] - boxCenter[i];
//  if (delta > 0.0) vectorB[i] = -1; // NOT
//  if (delta < 0.0) vectorB[i] = 1;  // NOT
//
  }
  mprintf("DEBUG: !VectorB = %2i %2i %2i\n",vectorB[0],vectorB[1],vectorB[2]);

  // A & !B
  ixyz[0]=0; ixyz[1]=0; ixyz[2]=0;
  for (i=0; i<3; i++) {
    if (vectorA[i] == vectorB[i]) ixyz[i] = vectorA[i];
    //ixyz[i] = vectorA[i] & vectorB[i];
  }
}
*/

/*
 * MinImageNonOrtho2()
 * Given two sets of coordinates and reciprocal space information based on
 * the current non-orthorhombic box, return the shortest imaged distance^2
 * between the coordinates.
 * The integer coefficients describing the closest reflection in reciprocal
 * space will be placed in ixyz.
 */
double MinImageNonOrtho2(double *Coord1, double *Coord2, double *box, bool origin, int *ixyz,
                         double *ucell, double *recip) {
  double min, f[3], f2[3];

  min = 100.0 * (box[0]*box[0]+box[1]*box[1]+box[2]*box[2]);

  //if (prnlev > 6) {
  //  fprintf(stdout, "ATOM      0  XXX A1      1     %7.3f %7.3f %7.3f\n",
  //          x1, y1, z1);
  //  fprintf(stdout, "ATOM      1  XXX A2      1     %7.3f %7.3f %7.3f\n",
  //          x2, y2, z2);
  //}

  f[0] = Coord1[0]*recip[0] + Coord1[1]*recip[1] + Coord1[2]*recip[2];
  f[1] = Coord1[0]*recip[3] + Coord1[1]*recip[4] + Coord1[2]*recip[5];
  f[2] = Coord1[0]*recip[6] + Coord1[1]*recip[7] + Coord1[2]*recip[8];

  f2[0] = Coord2[0]*recip[0] + Coord2[1]*recip[1] + Coord2[2]*recip[2];
  f2[1] = Coord2[0]*recip[3] + Coord2[1]*recip[4] + Coord2[2]*recip[5];
  f2[2] = Coord2[0]*recip[6] + Coord2[1]*recip[7] + Coord2[2]*recip[8];

  if (origin) {
    f[0] += 0.5;
    f[1] += 0.5;
    f[2] += 0.5;
    f2[0] += 0.5;
    f2[1] += 0.5;
    f2[2] += 0.5;
  }

  min = DIST2_ImageNonOrtho(f, f2, min, ixyz, ucell);

  return min;
}

/*
 * Frame::DIST2_ImageNonOrtho()
 * Given two coordinates and reciprocal space information based on 
 * the current non-orthorhombic box, return the shortest imaged distance^2 
 * between the coordinates.
 */
double DIST2_ImageNonOrtho(double *a1, double *a2, double *ucell, double *recip) { 
// double closest2
  double f[3], f2[3];
  int ixyz[3];

  f[0] = a2[0]*recip[0] + a2[1]*recip[1] + a2[2]*recip[2];
  f[1] = a2[0]*recip[3] + a2[1]*recip[4] + a2[2]*recip[5];
  f[2] = a2[0]*recip[6] + a2[1]*recip[7] + a2[2]*recip[8];

  f2[0] = a1[0]*recip[0] + a1[1]*recip[1] + a1[2]*recip[2];
  f2[1] = a1[0]*recip[3] + a1[1]*recip[4] + a1[2]*recip[5];
  f2[2] = a1[0]*recip[6] + a1[1]*recip[7] + a1[2]*recip[8];

  return DIST2_ImageNonOrtho(f, f2, -1.0, ixyz, ucell);
}

/*
 * DIST2_ImageNonOrtho()
 * Given two coordinate sets in reciprocal space, return the minimum imaged
 * distance^2 between them.
 * If minIn is > 0.0 it is considered a possible minimum distance.
 * The integer coefficients describing the closest reflection in reciprocal
 * space will be placed in ixyz.
 */
double DIST2_ImageNonOrtho(double *f, double *f2, double minIn, int *ixyz, double *ucell) { 
// double closest2
  double fx, fy, fz, f2x, f2y, f2z;
  double x,y,z,D,min;
  int ix,iy,iz;
  // DEBUG

  /*
   *  NON-ORTHORHOMBIC CASE: find shortest distance in periodic reference
   *  This is a brute force check requiring up to 26 distance evaluations.
   *  It has been adapted to be smarter by returning the first distance that
   *  is shorter than the minimum possible distance between images.
   */

  fx = f[0] - floor(f[0]);
  fy = f[1] - floor(f[1]);
  fz = f[2] - floor(f[2]); 
  
  f2x = f2[0] - floor(f2[0]);
  f2y = f2[1] - floor(f2[1]);
  f2z = f2[2] - floor(f2[2]);

  // Calc ix iy iz = 0 case
  x = (fx*ucell[0] + fy*ucell[3] + fz*ucell[6]) - (f2x*ucell[0] + f2y*ucell[3] + f2z*ucell[6]);
  y = (fx*ucell[1] + fy*ucell[4] + fz*ucell[7]) - (f2x*ucell[1] + f2y*ucell[4] + f2z*ucell[7]);
  z = (fx*ucell[2] + fy*ucell[5] + fz*ucell[8]) - (f2x*ucell[2] + f2y*ucell[5] + f2z*ucell[8]);
  // DEBUG
  //fprintf(stdout,"DEBUG: a2: fx  fy  fz  = %lf %lf %lf\n",fx,fy,fz);
  //fprintf(stdout,"DEBUG: a1: f2x f2y f2z = %lf %lf %lf\n",f2x,f2y,f2z);
  min = (x*x) + (y*y) + (z*z);

  if (minIn > 0.0 && minIn < min) min = minIn;

  ixyz[0] = 0;
  ixyz[1] = 0;
  ixyz[2] = 0;

  //if (closest2 != 0.0 && min < closest2) return (min);
//  this->ClosestImage(a1, a2, ixyz);
//  fprintf(stdout,"DEBUG: Predict  = %2i %2i %2i\n",ixyz[0],ixyz[1],ixyz[2]);

//  ix = ixyz[0];
//  iy = ixyz[1];
//  iz = ixyz[2];

  for (ix = -1; ix <= 1; ix++) {
    for (iy = -1; iy <= 1; iy++) {
      for (iz = -1; iz <= 1; iz++) {

        if (! (ix == 0 && iy == 0 && iz == 0) ) {
          x = ((fx+ix)*ucell[0] + (fy+iy)*ucell[3] + (fz+iz)*ucell[6]) - 
              (f2x*ucell[0] + f2y*ucell[3] + f2z*ucell[6]);
          y = ((fx+ix)*ucell[1] + (fy+iy)*ucell[4] + (fz+iz)*ucell[7]) - 
              (f2x*ucell[1] + f2y*ucell[4] + f2z*ucell[7]);
          z = ((fx+ix)*ucell[2] + (fy+iy)*ucell[5] + (fz+iz)*ucell[8]) - 
              (f2x*ucell[2] + f2y*ucell[5] + f2z*ucell[8]);
          D = (x*x) + (y*y) + (z*z);

          //if (debug > 3) 
          //  printf("DISTANCE + %2i*X %2i*Y %2i*Z unit cells is %8.3f\n", ix, iy, iz, sqrt(D));
          
          if (D < min) {
            min = D;
            ixyz[0] = ix;
            ixyz[1] = iy;
            ixyz[2] = iz;
            // DEBUG
            //if (closest2 != 0.0 && min < closest2)
            //  return(min);
          }
        }

      }
    }
  }

  //D = sqrt(min);
//  fprintf(stdout,"DEBUG: MinDist  = %2i %2i %2i = %8.3f\n", ixmin, iymin, izmin, D);
//  printf("---------------------------------------------------------------\n");
  return(min);
}

/*
 * Frame::DIST2_ImageOrtho()
 * Return the minimum orthorhombic imaged distance^2 between coordinates a1 
 * and a2.
 */
double DIST2_ImageOrtho(double *a1, double *a2, double *box) {
  double x,y,z,D;

  x = a1[0] - a2[0];
  y = a1[1] - a2[1];
  z = a1[2] - a2[2];

  // Get rid of sign info
  if (x<0) x=-x;
  if (y<0) y=-y;
  if (z<0) z=-z;

  // Get rid of multiples of box lengths 
  while (x > box[0]) x = x - box[0];
  while (y > box[1]) y = y - box[1];
  while (z > box[2]) z = z - box[2];

  // Find shortest distance in periodic reference
  D = box[0] - x;
  if (D < x) x = D;
  D = box[1] - y;
  if (D < y) y = D;  
  D = box[0] - z;
  if (D < z) z = D;

  x = x * x;
  y = y * y;
  z = z * z;
 
  //D = sqrt(x + y + z);
  D = x + y + z;

  return D;
}

/*
 * Frame::DIST2_NoImage()
 * Return distance^2 between coordinates in a1 and a2.
 */
double DIST2_NoImage(double *a1, double *a2) {
  double x,y,z,D;

  x = a1[0] - a2[0];
  y = a1[1] - a2[1];
  z = a1[2] - a2[2];

  x=x*x;
  y=y*y;
  z=z*z;

  //D=sqrt(x + y + z);
  D = x + y + z;

  //fprintf(stdout,"Mask1=%8.3lf %8.3lf %8.3lf Mask2=%8.3lf %8.3lf %8.3lf D=%8.3lf\n",
  //        a1[0],a1[1],a1[2],a2[0],a2[1],a2[2],D);

  return D;
}


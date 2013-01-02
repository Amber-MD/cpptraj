#include <cmath> //floor
#include "DistRoutines.h"

// Frame::DIST2_ImageNonOrtho()
/** Given two coordinates and reciprocal space information based on 
  * the current non-orthorhombic box, return the shortest imaged distance^2 
  * between the coordinates.
  */
double DIST2_ImageNonOrtho(Vec3 const& a1, Vec3 const& a2, 
                           Matrix_3x3 const& ucell, Matrix_3x3 const& recip) 
{ 
  int ixyz[3];
  Vec3 f = recip * a2;
  Vec3 f2 = recip * a1;

  return DIST2_ImageNonOrthoRecip(f, f2, -1.0, ixyz, ucell);
}

// DIST2_ImageNonOrthoRecip()
/** Given two coordinate sets in reciprocal space, return the minimum imaged
  * distance^2 between them.
  * If minIn is > 0.0 it is considered a possible minimum distance.
  * The integer coefficients describing the closest reflection in reciprocal
  * space will be placed in ixyz.
  */
double DIST2_ImageNonOrthoRecip(Vec3 const& f, Vec3 const& f2, double minIn, 
                                int* ixyz, Matrix_3x3 const& ucell) 
{ 
  //double closest2
  double fx, fy, fz, f2x, f2y, f2z, X_factor, Y_factor, Z_factor;
  double fxm1, fxp1, fym1, fyp1, fzm1, fzp1;
  double x,y,z,D,min;

  /*  NON-ORTHORHOMBIC CASE: find shortest distance in periodic reference
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

  // Precompute some factors
  X_factor = (f2x*ucell[0] + f2y*ucell[3] + f2z*ucell[6]);
  Y_factor = (f2x*ucell[1] + f2y*ucell[4] + f2z*ucell[7]);
  Z_factor = (f2x*ucell[2] + f2y*ucell[5] + f2z*ucell[8]);

  fxm1 = fx - 1; fxp1 = fx + 1;
  fym1 = fy - 1; fyp1 = fy + 1;
  fzm1 = fz - 1; fzp1 = fz + 1;

  double fxm1u0 = fxm1 * ucell[0];
  double fxu0   = fx   * ucell[0];
  double fxp1u0 = fxp1 * ucell[0];
  double fxm1u1 = fxm1 * ucell[1];
  double fxu1   = fx   * ucell[1];
  double fxp1u1 = fxp1 * ucell[1];
  double fxm1u2 = fxm1 * ucell[2];
  double fxu2   = fx   * ucell[2];
  double fxp1u2 = fxp1 * ucell[2];

  double fym1u3 = fym1 * ucell[3];
  double fyu3   = fy   * ucell[3];
  double fyp1u3 = fyp1 * ucell[3];
  double fym1u4 = fym1 * ucell[4];
  double fyu4   = fy   * ucell[4];
  double fyp1u4 = fyp1 * ucell[4];
  double fym1u5 = fym1 * ucell[5];
  double fyu5   = fy   * ucell[5];
  double fyp1u5 = fyp1 * ucell[5];

  double fzm1u6 = fzm1 * ucell[6];
  double fzu6   = fz   * ucell[6];
  double fzp1u6 = fzp1 * ucell[6];
  double fzm1u7 = fzm1 * ucell[7];
  double fzu7   = fz   * ucell[7];
  double fzp1u7 = fzp1 * ucell[7];
  double fzm1u8 = fzm1 * ucell[8];
  double fzu8   = fz   * ucell[8];
  double fzp1u8 = fzp1 * ucell[8];

  // Calc ix iy iz = 0 case
  x = (fxu0 + fyu3 + fzu6) - X_factor;
  y = (fxu1 + fyu4 + fzu7) - Y_factor;
  z = (fxu2 + fyu5 + fzu8) - Z_factor;
  // DEBUG
  //fprintf(stdout,"DEBUG: a2: fx  fy  fz  = %lf %lf %lf\n",fx,fy,fz);
  //fprintf(stdout,"DEBUG: a1: f2x f2y f2z = %lf %lf %lf\n",f2x,f2y,f2z);
  min = (x*x) + (y*y) + (z*z);

  if (minIn > 0.0 && minIn < min) min = minIn;

  ixyz[0] = 0;
  ixyz[1] = 0;
  ixyz[2] = 0;

  // -1 -1 -1
  x = (fxm1u0 + fym1u3 + fzm1u6) - X_factor;
  y = (fxm1u1 + fym1u4 + fzm1u7) - Y_factor;
  z = (fxm1u2 + fym1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] = -1; ixyz[2] = -1; }
  // -1 -1  0
  x = (fxm1u0 + fym1u3 + fzu6  ) - X_factor;
  y = (fxm1u1 + fym1u4 + fzu7  ) - Y_factor;
  z = (fxm1u2 + fym1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] = -1; ixyz[2] =  0; }
  // -1 -1 +1
  x = (fxm1u0 + fym1u3 + fzp1u6) - X_factor;
  y = (fxm1u1 + fym1u4 + fzp1u7) - Y_factor;
  z = (fxm1u2 + fym1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] = -1; ixyz[2] =  1; }
  // -1  0 -1
  x = (fxm1u0 + fyu3   + fzm1u6) - X_factor;
  y = (fxm1u1 + fyu4   + fzm1u7) - Y_factor;
  z = (fxm1u2 + fyu5   + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  0; ixyz[2] = -1; }
  // -1  0  0
  x = (fxm1u0 + fyu3   + fzu6  ) - X_factor;
  y = (fxm1u1 + fyu4   + fzu7  ) - Y_factor;
  z = (fxm1u2 + fyu5   + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  0; ixyz[2] =  0; }
  // -1  0 +1
  x = (fxm1u0 + fyu3   + fzp1u6) - X_factor;
  y = (fxm1u1 + fyu4   + fzp1u7) - Y_factor;
  z = (fxm1u2 + fyu5   + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  0; ixyz[2] =  1; }
  // -1 +1 -1
  x = (fxm1u0 + fyp1u3 + fzm1u6) - X_factor;
  y = (fxm1u1 + fyp1u4 + fzm1u7) - Y_factor;
  z = (fxm1u2 + fyp1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  1; ixyz[2] = -1; }
  // -1 +1  0
  x = (fxm1u0 + fyp1u3 + fzu6  ) - X_factor;
  y = (fxm1u1 + fyp1u4 + fzu7  ) - Y_factor;
  z = (fxm1u2 + fyp1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  1; ixyz[2] =  0; }
  // -1 +1 +1
  x = (fxm1u0 + fyp1u3 + fzp1u6) - X_factor;
  y = (fxm1u1 + fyp1u4 + fzp1u7) - Y_factor;
  z = (fxm1u2 + fyp1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  1; ixyz[2] =  1; }

  //  0 -1 -1
  x = (fxu0   + fym1u3 + fzm1u6) - X_factor;
  y = (fxu1   + fym1u4 + fzm1u7) - Y_factor;
  z = (fxu2   + fym1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] = -1; ixyz[2] = -1; }
  //  0 -1  0
  x = (fxu0   + fym1u3 + fzu6  ) - X_factor;
  y = (fxu1   + fym1u4 + fzu7  ) - Y_factor;
  z = (fxu2   + fym1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] = -1; ixyz[2] =  0; }
  //  0 -1 +1
  x = (fxu0   + fym1u3 + fzp1u6) - X_factor;
  y = (fxu1   + fym1u4 + fzp1u7) - Y_factor;
  z = (fxu2   + fym1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] = -1; ixyz[2] =  1; }
  //  0  0 -1
  x = (fxu0   + fyu3   + fzm1u6) - X_factor;
  y = (fxu1   + fyu4   + fzm1u7) - Y_factor;
  z = (fxu2   + fyu5   + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  0; ixyz[2] = -1; }
  //  0  0  0
  //  0  0 +1
  x = (fxu0   + fyu3   + fzp1u6) - X_factor;
  y = (fxu1   + fyu4   + fzp1u7) - Y_factor;
  z = (fxu2   + fyu5   + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  0; ixyz[2] =  1; }
  //  0 +1 -1
  x = (fxu0   + fyp1u3 + fzm1u6) - X_factor;
  y = (fxu1   + fyp1u4 + fzm1u7) - Y_factor;
  z = (fxu2   + fyp1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  1; ixyz[2] = -1; }
  //  0 +1  0
  x = (fxu0   + fyp1u3 + fzu6  ) - X_factor;
  y = (fxu1   + fyp1u4 + fzu7  ) - Y_factor;
  z = (fxu2   + fyp1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  1; ixyz[2] =  0; }
  //  0 +1 +1
  x = (fxu0   + fyp1u3 + fzp1u6) - X_factor;
  y = (fxu1   + fyp1u4 + fzp1u7) - Y_factor;
  z = (fxu2   + fyp1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  1; ixyz[2] =  1; }

  // +1 -1 -1
  x = (fxp1u0 + fym1u3 + fzm1u6) - X_factor;
  y = (fxp1u1 + fym1u4 + fzm1u7) - Y_factor;
  z = (fxp1u2 + fym1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] = -1; ixyz[2] = -1; }
  // +1 -1  0
  x = (fxp1u0 + fym1u3 + fzu6  ) - X_factor;
  y = (fxp1u1 + fym1u4 + fzu7  ) - Y_factor;
  z = (fxp1u2 + fym1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] = -1; ixyz[2] =  0; }
  // +1 -1 +1
  x = (fxp1u0 + fym1u3 + fzp1u6) - X_factor;
  y = (fxp1u1 + fym1u4 + fzp1u7) - Y_factor;
  z = (fxp1u2 + fym1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] = -1; ixyz[2] =  1; }
  // +1  0 -1
  x = (fxp1u0 + fyu3   + fzm1u6) - X_factor;
  y = (fxp1u1 + fyu4   + fzm1u7) - Y_factor;
  z = (fxp1u2 + fyu5   + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  0; ixyz[2] = -1; }
  // +1  0  0
  x = (fxp1u0 + fyu3   + fzu6  ) - X_factor;
  y = (fxp1u1 + fyu4   + fzu7  ) - Y_factor;
  z = (fxp1u2 + fyu5   + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  0; ixyz[2] =  0; }
  // +1  0 +1
  x = (fxp1u0 + fyu3   + fzp1u6) - X_factor;
  y = (fxp1u1 + fyu4   + fzp1u7) - Y_factor;
  z = (fxp1u2 + fyu5   + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  0; ixyz[2] =  1; }
  // +1 +1 -1
  x = (fxp1u0 + fyp1u3 + fzm1u6) - X_factor;
  y = (fxp1u1 + fyp1u4 + fzm1u7) - Y_factor;
  z = (fxp1u2 + fyp1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  1; ixyz[2] = -1; }
  // +1 +1  0
  x = (fxp1u0 + fyp1u3 + fzu6  ) - X_factor;
  y = (fxp1u1 + fyp1u4 + fzu7  ) - Y_factor;
  z = (fxp1u2 + fyp1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  1; ixyz[2] =  0; }
  // +1 +1 +1
  x = (fxp1u0 + fyp1u3 + fzp1u6) - X_factor;
  y = (fxp1u1 + fyp1u4 + fzp1u7) - Y_factor;
  z = (fxp1u2 + fyp1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  1; ixyz[2] =  1; }

  //if (closest2 != 0.0 && min < closest2) return (min);
//  this->ClosestImage(a1, a2, ixyz);
//  fprintf(stdout,"DEBUG: Predict  = %2i %2i %2i\n",ixyz[0],ixyz[1],ixyz[2]);

//  ix = ixyz[0];
//  iy = ixyz[1];
//  iz = ixyz[2];

//D = sqrt(min);
//  fprintf(stdout,"DEBUG: MinDist  = %2i %2i %2i = %8.3f\n", ixmin, iymin, izmin, D);
//  printf("---------------------------------------------------------------\n");
  return(min);
}

// Frame::DIST2_ImageOrtho()
/** Return the minimum orthorhombic imaged distance^2 between coordinates a1 
  * and a2.
  */
double DIST2_ImageOrtho(Vec3 const& a1, Vec3 const& a2, Box const& box) {
  // If box lengths are zero no imaging possible
  if (box[0]==0.0 && box[1]==0.0 && box[2]==0.0) return -1.0;
  double x = a1[0] - a2[0];
  double y = a1[1] - a2[1];
  double z = a1[2] - a2[2];
  // Get rid of sign info
  if (x<0) x=-x;
  if (y<0) y=-y;
  if (z<0) z=-z;
  // Get rid of multiples of box lengths 
  while (x > box[0]) x = x - box[0];
  while (y > box[1]) y = y - box[1];
  while (z > box[2]) z = z - box[2];
  // Find shortest distance in periodic reference
  double D = box[0] - x;
  if (D < x) x = D;
  D = box[1] - y;
  if (D < y) y = D;  
  D = box[2] - z;
  if (D < z) z = D;

  return (x*x + y*y + z*z);
}

// Frame::DIST2_NoImage()
/** Return distance^2 between coordinates in a1 and a2.
  */
double DIST2_NoImage(const double* a1, const double* a2) {
  double x = a1[0] - a2[0];
  double y = a1[1] - a2[1];
  double z = a1[2] - a2[2];
  //double D = x*x + y*y + z*z;
  //fprintf(stdout,"Mask1=%8.3f %8.3f %8.3f Mask2=%8.3f %8.3f %8.3f D=%8.3f\n",
  //        a1[0],a1[1],a1[2],a2[0],a2[1],a2[2],D);
  return (x*x + y*y + z*z);
}

double DIST2_NoImage( Vec3 const& a1, Vec3 const& a2 ) {
  Vec3 vec = a1 - a2;
  return vec.Magnitude2();
}

double DIST2(const double* a1, const double* a2, ImagingType itype,
             Box const& box, Matrix_3x3 const& ucell, Matrix_3x3 const& recip)
{
  if (itype==NOIMAGE) 
    return DIST2_NoImage( a1, a2 );
  else if (itype==ORTHO) 
    return DIST2_ImageOrtho( a1, a2, box );
  else // NONORTHO
    return DIST2_ImageNonOrtho( a1, a2, ucell, recip );
}

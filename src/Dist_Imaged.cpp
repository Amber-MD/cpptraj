#include "Dist_Imaged.h"
//#incl ude "CpptrajStdio.h" // DEBUG
#include "Matrix_3x3.h"
#include "Vec3.h"
#include <limits> // numeric_limits<double>::max()
//#incl ude <cmath> // DEBUG

/// Bring each coordinate back into the main unit cell (0 to 1)
static inline void wrap_frac(Vec3& frac) {
  for (int idx = 0; idx != 3; idx++) {
    while (frac[idx] < 0.0)
      frac[idx] += 1.0;
    while (frac[idx] > 1.0)
      frac[idx] -= 1.0;
  }
}

/** Calculate minimum imaged distance^2 between two points in Cartesian space.
  * Also set the indices of the closest image.
  */
double Cpptraj::Dist2_Imaged(Vec3 const& xyz1, Vec3 const& xyz2,
                             Matrix_3x3 const& ucell, Matrix_3x3 const& fcell,
                             int* ixyz)
{
  // Wrap each coordinate in fractional coordinates back to the main unit cell
  // (so each frac coord is between 0 and 1).
  Vec3 f1 = fcell * xyz1;

  Vec3 f2 = fcell * xyz2;

  wrap_frac(f1);
  wrap_frac(f2);

//  // DEBUG: f1 to f2 in fractional space
//  Vec3 fdelta = f2 - f1;
//  Vec3 fdelta_abs = fdelta;
//  for (int idx = 0; idx != 3; idx++) // DEBUG
//    if (fdelta_abs[idx] < 0) fdelta_abs[idx] = -fdelta_abs[idx]; // DEBUG

  // Calculate wrapped f2 back in Cartesian space
  Vec3 c2 = ucell.TransposeMult( f2 );

  // Brute force case. Check each image of f1 to determine which distance is the smallest.
  double minD2 = std::numeric_limits<double>::max();
  ixyz[0] = 0;
  ixyz[1] = 0;
  ixyz[2] = 0;
  for (int ix = -1; ix != 2; ix++) {
    for (int iy = -1; iy != 2; iy++) {
      for (int iz = -1; iz != 2; iz++) {
        // Image in frac space
        Vec3 f1img = f1 + Vec3(ix, iy, iz);
        // Image in Cartesian space
        Vec3 c1img = ucell.TransposeMult( f1img );
        // Calculate distance in Cartesian space
        double dx = (c1img[0] - c2[0]);
        double dy = (c1img[1] - c2[1]);
        double dz = (c1img[2] - c2[2]);
        double D2 = (dx*dx) + (dy*dy) + (dz*dz);
        if (D2 < minD2) {
          minD2 = D2;
          ixyz[0] = ix;
          ixyz[1] = iy;
          ixyz[2] = iz;
        }
      } // END loop over Z image indices
    } // END loop over Y image indices
  } // END loop over X image indices

  // DEBUG
//  for (int idx = 0; idx != 3; idx++) {
//    if (ixyz[idx] != 0) {
//      int fdelta_sign = 0;
//      if (fdelta[idx] > 0)
//        fdelta_sign = 1;
//      else if (fdelta[idx] < 0)
//        fdelta_sign = -1;
//      if (fdelta_sign != ixyz[idx]) {
//    //if (ixyz[idx] != 0 && fdelta_abs[idx] < 0.5)
//        mprintf("DEBUG: ixyz %2i %2i %2i  fdelta= %12.4f %12.4f %12.4f\n",
//                ixyz[0], ixyz[1], ixyz[2], fdelta[0], fdelta[1], fdelta[2]);
//      }
//    }
//  } 

  return minD2;
}

/** Calculate minimum imaged distance^2 between two points in Fractional space.
  * Assumes fractional coordinates have already been wrapped into the 
  * primary unit cell. Also set the indices of the closest image of
  * f1 to f2.
  * TODO use in place of DIST2_ImageNonOrthoRecip()
  */
double Cpptraj::Dist2_Imaged_Frac(Vec3 const& f1, Vec3 const& f2,
                                  Matrix_3x3 const& ucell, Matrix_3x3 const& fcell,
                                  int* ixyz)
{
  // Calculate f2 to f1 in fractional space
  Vec3 fdelta = f2 - f1;
  //fdelta.Print("fdelta");
  // Set the min and max cell indices we need to look for based
  // on the f1->f2 vector.
  int imin[3];
  int imax[3];
  for (int idx = 0; idx != 3; idx++) {
    if (fdelta[idx] < 0) {
      imin[idx] = -1;
      imax[idx] = 1;
    } else {
      imin[idx] = 0;
      imax[idx] = 2;
    }
  }


  // Calculate wrapped f2 back in Cartesian space
  Vec3 c2 = ucell.TransposeMult( f2 );
/*
  // Set the min and max cell indices we need to look for based on 
  // where f2 is located.
  int imin[3];
  int imax[3];
  for (int idx = 0; idx != 3; idx++) {
    if (f2[idx] < 0.5) {
      // Need to check in the negative direction
      imin[idx] = -1;
      imax[idx] = 1;
    } else if (f2[idx] > 0.5) {
      // Need to check in the positive direction
      imin[idx] = 0;
      imax[idx] = 2;
    } else {
      // Exact center - check everywhere
      imin[idx] = -1;
      imax[idx] = 2;
    }
  }*/

  // Searching only images we are close to.
  // Check each image of f1 to determine which distance is the smallest.
  double minD2 = std::numeric_limits<double>::max();
  ixyz[0] = 0;
  ixyz[1] = 0;
  ixyz[2] = 0;
  for (int ix = imin[0]; ix != imax[0]; ix++) {
    for (int iy = imin[1]; iy != imax[1]; iy++) {
      for (int iz = imin[2]; iz != imax[2]; iz++) {
        // Image in frac space
        Vec3 f1img = f1 + Vec3(ix, iy, iz);
        // Image in Cartesian space
        Vec3 c1img = ucell.TransposeMult( f1img );
        // Calculate distance in Cartesian space
        double dx = (c1img[0] - c2[0]);
        double dy = (c1img[1] - c2[1]);
        double dz = (c1img[2] - c2[2]);
        double D2 = (dx*dx) + (dy*dy) + (dz*dz);
        if (D2 < minD2) {
          minD2 = D2;
          ixyz[0] = ix;
          ixyz[1] = iy;
          ixyz[2] = iz;
        }
      } // END loop over Z image indices
    } // END loop over Y image indices
  } // END loop over X image indices
  //mprintf("minD %12.4f\n", sqrt(minD2));
  return minD2;
}

/** Calculate minimum imaged distance^2 between two points in Cartesian space.
  * Also set the indices of the closest image of xyz1 to xyz2.
  * TODO use in place of DIST2_ImageNonOrtho()
  */
double Cpptraj::Dist2_Imaged_Cart(Vec3 const& xyz1, Vec3 const& xyz2,
                                  Matrix_3x3 const& ucell, Matrix_3x3 const& fcell,
                                  int* ixyz)
{
  // Wrap each coordinate in fractional coordinates back to the main unit cell
  // (so each frac coord is between 0 and 1).
  Vec3 f1 = fcell * xyz1;

  Vec3 f2 = fcell * xyz2;

  wrap_frac(f1);
  wrap_frac(f2);

  return Dist2_Imaged_Frac(f1, f2, ucell, fcell, ixyz);
}

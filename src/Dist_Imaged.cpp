#include "Dist_Imaged.h"
#include "Matrix_3x3.h"
#include "Vec3.h"
#include <limits>

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
  return minD2;
}

#ifndef INC_UNWRAP_H
#define INC_UNWRAP_H
#include "Matrix_3x3.h"
#include "Vec3.h"
namespace Cpptraj {
namespace Unwrap {
  /// \return True if absolute value of any component of dxyz is greater than limxyz
  template <typename T> bool vGreaterThan(T const* dxyz, T const* limxyz) {
    return ( fabs(dxyz[0]) > limxyz[0] ||
             fabs(dxyz[1]) > limxyz[1] ||
             fabs(dxyz[2]) > limxyz[2] );
  }

  /// \return Vector required for unwrapping in orthogonal box
  template <typename T> Vec3 UnwrapVec_Ortho(Vec3 const& vtgt, Vec3 const& vref, Vec3 const& boxVec, Vec3 const& limxyz)
  {
    T dxyz[3];
    dxyz[0] = (vtgt[0] - vref[0]);
    dxyz[1] = (vtgt[1] - vref[1]);
    dxyz[2] = (vtgt[2] - vref[2]);
    if (vGreaterThan<T>(dxyz, limxyz.Dptr()))
      return Vec3(
        -floor( dxyz[0] / boxVec[0] + 0.5 ) * boxVec[0],
        -floor( dxyz[1] / boxVec[1] + 0.5 ) * boxVec[1],
        -floor( dxyz[2] / boxVec[2] + 0.5 ) * boxVec[2]
      );
    else
      return Vec3(0.0);
  }

  /// \return Vector required for unwrapping in non-orthogonal box
  template <typename T> Vec3 UnwrapVec_Nonortho(Vec3 const& vtgt, Vec3 const& vref, Matrix_3x3 const& ucell, Matrix_3x3 const& frac, Vec3 const& limxyz)
  {
    Vec3 boxTrans(0.0);
    // Calculate original distance from the ref (previous) position. 
    Vec3 vd = vtgt - vref; // dx dy dz
    if (!vGreaterThan<T>(vd.Dptr(), limxyz.Dptr()))
      return boxTrans;
    T minDistanceSquare = vd.Magnitude2();
    // Reciprocal coordinates
    vd = frac * vd ; // recip * dxyz
    T cx = floor(vd[0]);
    T cy = floor(vd[1]);
    T cz = floor(vd[2]);
    // Loop over all possible translations 
    for (int ix = -1; ix < 2; ++ix) {
      for (int iy = -1; iy < 2; ++iy) {
        for (int iz = -1; iz < 2; ++iz) {
          // Calculate the translation.
          Vec3 vcc = ucell.TransposeMult( Vec3( cx+(T)ix, 
                                                cy+(T)iy, 
                                                cz+(T)iz ) ); // ucell^T * ccxyz
          // Calc. the potential new coordinate for tgt
          Vec3 vnew = vtgt - vcc; 
          // Calc. the new distance from the ref (previous) position
          Vec3 vr = vref - vnew; 
          T distanceSquare = vr.Magnitude2();
          // If the orig. distance is greater than the new distance, unwrap. 
          if ( minDistanceSquare > distanceSquare ) {
              minDistanceSquare = distanceSquare;
              boxTrans = vcc;
          }
        }
      }
    }
    // Translate tgt atoms
    boxTrans.Neg();
    return boxTrans;
  }

} // END namespace Unwrap
} // END namespace Cpptraj
#endif

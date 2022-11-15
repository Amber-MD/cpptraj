#ifndef INC_ENERGY_ENE_BOND_H
#define INC_ENERGY_ENE_BOND_H
#include "Kernel_Harmonic.h"
namespace Cpptraj {
namespace Energy {
/// \return Energy of the bond between xyz1 and xyz2
template <typename T>
T Ene_Bond(T const* xyz1, T const* xyz2, T const& req, T const& rk)
{
  T x = xyz1[0] - xyz2[0];
  T y = xyz1[1] - xyz2[1];
  T z = xyz1[2] - xyz2[2];
  T r2 = (x*x + y*y + z*z);
  T r = sqrt(r2);
  T ene = Kernel_Harmonic<T>( r, rk, req );
  return ene;
}
}
}
#endif

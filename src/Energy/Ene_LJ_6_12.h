#ifndef INC_ENERGY_ENE_LJ_6_12_H
#define INC_ENERGY_ENE_LJ_6_12_H
namespace Cpptraj {
namespace Energy {
/// \return LJ 6-12 energy
template <typename T>
T Ene_LJ_6_12(T const& rij2, T const& LJA, T const& LJB)
{
  T r2    = 1.0 / rij2;
  T r6    = r2  * r2 * r2;
  T r12   = r6  * r6;
  T f12   = LJA * r12;  // A/r^12
  T f6    = LJB * r6;   // B/r^6
  T e_vdw = f12 - f6;   // (A/r^12)-(B/r^6)
  return e_vdw;
}
}
}
#endif

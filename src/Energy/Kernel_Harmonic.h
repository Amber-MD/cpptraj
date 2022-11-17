#ifndef INC_ENERGYKERNEL_HARMONIC_H
#define INC_ENERGYKERNEL_HARMONIC_H
namespace Cpptraj {
namespace Energy {
/// \return Hooke's law type energy
template <typename T> T Kernel_Harmonic(T r, T Rk, T Req) {
  T rdiff = r - Req;
  T ene = Rk * (rdiff * rdiff);
  return ene;
}
}
}
#endif

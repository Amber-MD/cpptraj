#ifndef INC_ENERGYKERNEL_HARMONIC_H
#define INC_ENERGYKERNEL_HARMONIC_H
template <typename T> T EnergyKernel_Harmonic(T r, T Rk, T Req) {
  T rdiff = r - Req;
  T ene = Rk * (rdiff * rdiff);
  return ene;
};
#endif

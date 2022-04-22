#ifndef INC_ENERGYKERNEL_HARMONICBOND_H
#define INC_ENERGYKERNEL_HARMONICBOND_H
template <typename T> T EnergyKernel_HarmonicBond(T r, T Rk, T Req) {
  T rdiff = r - Req;
  T ene = Rk * (rdiff * rdiff);
  return ene;
};
#endif

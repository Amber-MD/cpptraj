#ifndef INC_ENERGYKERNEL_FOURIER_H
#define INC_ENERGYKERNEL_FOURIER_H
template <typename T> T EnergyKernel_Fourier(T phi, T Pk, T Pn, T Phase) {
  T ene = Pk * (1 + cos(Pn * phi - Phase));
  return ene;
}
#endif

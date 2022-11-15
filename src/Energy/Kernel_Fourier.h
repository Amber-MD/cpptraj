#ifndef INC_ENERGYKERNEL_FOURIER_H
#define INC_ENERGYKERNEL_FOURIER_H
namespace Cpptraj {
namespace Energy {
template <typename T> T Kernel_Fourier(T phi, T Pk, T Pn, T Phase) {
  T ene = Pk * (1 + cos(Pn * phi - Phase));
  return ene;
}
}
}
#endif

#ifndef INC_KERNEL_EWALDADJUST_H
#define INC_KERNEL_EWALDADJUST_H
namespace Cpptraj {
namespace Energy {
/// Used to adjust Ewald energy for excluded atoms
template <typename T> T Kernel_EwaldAdjust(T const& q0, T const& q1, T const& rij, T const& erfcval)
{
  T d0 = (erfcval - 1.0) / rij;
  return (q0 * q1 * d0);
}
}
}
#endif

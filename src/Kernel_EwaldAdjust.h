#ifndef INC_KERNEL_EWALDADJUST_H
#define INC_KERNEL_EWALDADJUST_H
// NOTE: ERFC is not const so that timing data can be obtained
template <typename T> T Kernel_EwaldAdjust(T const& q0, T const& q1, T const& rij,
                                           T const& ew_coeff, ErfcFxn& ERFC)
{
  T erfc = ERFC.ERFC(ew_coeff * rij);
  T d0 = (erfc - 1.0) / rij;
  return (q0 * q1 * d0);
}
#endif

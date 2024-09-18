#ifndef INC_LJSWITCH_H
#define INC_LJSWITCH_H
namespace Cpptraj {
namespace Energy {
/** Switching function for Lennard-Jones. */
template <typename T> T LJswitch(T rij2, T cut2_0, T cut2_1)
{
  if (rij2 <= cut2_0)
    return 1.0;
  else if (rij2 > cut2_1)
    return 0.0;
  else {
    T xoff_m_x = cut2_1 - rij2;
    T fac = 1.0 / (cut2_1 - cut2_0);
    return (xoff_m_x*xoff_m_x) * (cut2_1 + 2.0*rij2 - 3.0*cut2_0) * (fac*fac*fac);
  }
}
}
}
#endif

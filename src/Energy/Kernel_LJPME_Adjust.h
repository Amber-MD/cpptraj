#ifndef INC_KERNEL_LJPME_ADJUST_H
#define INC_KERNEL_LJPME_ADJUST_H
namespace Cpptraj {
namespace Energy {
/// Used to adjust LJPME energy for excluded atoms
template <typename T> T Kernel_LJPME_Adjust(T const& rij2, T const& lj_ewcoeff, T const& Cij)
{
  // LJ PME direct space exclusion correction
  // NOTE: Assuming excluded pair is within cutoff
  T kr2 = lj_ewcoeff * lj_ewcoeff * rij2;
  T kr4 = kr2 * kr2;
  //double kr6 = kr2 * kr4;
  T expterm = exp(-kr2);
  T r4 = rij2 * rij2;
  T r6 = rij2 * r4;
  //T Cij = EW_.CalcCij(atom0.Idx(), atom1.Idx()); //Cparam_[it0->Idx()] * Cparam_[it1->Idx()];
  return (1.0 - (1.0 +  kr2 + kr4/2.0)*expterm) / r6 * Cij;
}
}
}
#endif

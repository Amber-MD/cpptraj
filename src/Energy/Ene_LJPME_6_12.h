#ifndef INC_ENERGY_ENE_LJPME_6_12_H
#define INC_ENERGY_ENE_LJPME_6_12_H
namespace Cpptraj {
namespace Energy {
/// \return LJ 6-12 energy
template <typename T>
void Ene_LJPME_6_12(T& e_vdw, T& e_pmevdw,
                    T const& rij2, T const& LJA, T const& LJB, T const& lj_ewcoeff, T const& Cij)
{
  T r2    = 1.0 / rij2;
  T r6    = r2  * r2 * r2;
  T r12   = r6  * r6;
  T f12   = LJA * r12;  // A/r^12
  T f6    = LJB * r6;   // B/r^6
    e_vdw = f12 - f6;   // (A/r^12)-(B/r^6)
  // LJ PME direct space correction
  T kr2 = lj_ewcoeff * lj_ewcoeff * rij2;
  T kr4 = kr2 * kr2;
  //double kr6 = kr2 * kr4;
  T expterm = exp(-kr2);
  //T Cij = EW_.CalcCij(atom0.Idx(), atom1.Idx()); //Cparam_[it0->Idx()] * Cparam_[it1->Idx()];
  //Eljpme_correction_ += (1.0 - (1.0 +  kr2 + kr4/2.0)*expterm) * r6 * vswitch * Cij;
  e_pmevdw = (1.0 - (1.0 +  kr2 + kr4/2.0)*expterm) * r6 * Cij;
}
}
}
#endif

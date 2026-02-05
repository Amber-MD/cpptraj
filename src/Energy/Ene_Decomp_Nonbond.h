#ifndef INC_ENERGY_ENE_DECOMP_NONBOND_H
#define INC_ENERGY_ENE_DECOMP_NONBOND_H
#include "Ene_LJ_6_12.h"
namespace Cpptraj {
namespace Energy {
/// Calculate simple decomposed nonbond energy with Coulomb and LJ 6-12
/** NOTE: To make this file more lightweight, there are no includes for
  *       Frame, CharMask, etc. It is assumed that before this file is
  *       included there will be at least '#include "Topology.h" and
  *       'include <cmath>' (for the sqrt).
  * Comple with -DCPPTRAJ_DEBUG_ENEDECOMP for more details on the individual
  * contributions.
  */
template <typename T>
void Ene_Decomp_Nonbond(Frame const& fIn, Topology const& tIn, CharMask const& selectedAtoms,
                        ExclusionArray const& Excluded, T const& QFAC,
                        T& Eelec, T& Evdw,
                        std::vector<double>& atom_elec,
                        std::vector<double>& atom_vdw)
{
  Eelec = 0;
  Evdw = 0;
  for (int atom1 = 0; atom1 < tIn.Natom(); atom1++) {
    bool atom1_is_selected = selectedAtoms.AtomInCharMask( atom1 );
    const double* crd1 = fIn.XYZ( atom1 );
    ExclusionArray::ExListType const& excludedAtoms = Excluded[atom1];
    for (int atom2 = atom1 + 1; atom2 < tIn.Natom(); atom2++) {
      bool atom2_is_selected = selectedAtoms.AtomInCharMask( atom2 );
      if (atom1_is_selected || atom2_is_selected) {
        ExclusionArray::ExListType::const_iterator it = excludedAtoms.find( atom2 );
        if (it == excludedAtoms.end()) {
          // Either atom1 or atom2 is selected and the interaction is not excluded.
          // TODO image distances?
          const double* crd2 = fIn.XYZ( atom2 );
          T dx = crd1[0] - crd2[0];
          T dy = crd1[1] - crd2[1];
          T dz = crd1[2] - crd2[2];
          T rij2 = (dx*dx) + (dy*dy) + (dz*dz);
          // VDW
          NonbondType const& LJ = tIn.GetLJparam(atom1, atom2);
          T e_vdw = Ene_LJ_6_12<T>( rij2, LJ.A(), LJ.B() );
#         ifdef CPPTRAJ_DEBUG_ENEDECOMP
          mprintf("DEBUG: VDW %f\n", e_vdw);
#         endif
          Evdw += e_vdw;
          T ene_half = e_vdw * 0.5;
          if (atom1_is_selected) atom_vdw[atom1] += ene_half;
          if (atom2_is_selected) atom_vdw[atom2] += ene_half;
          // Coulomb energy
          T rij = sqrt(rij2);
          T qiqj = QFAC * tIn[atom1].Charge() * tIn[atom2].Charge();
          T e_elec = qiqj / rij;
#         ifdef CPPTRAJ_DEBUG_ENEDECOMP
          mprintf("DEBUG: ELE %f\n", e_elec);
#         endif
          Eelec += e_elec;
          ene_half = e_elec * 0.5;
          if (atom1_is_selected) atom_elec[atom1] += ene_half;
          if (atom2_is_selected) atom_elec[atom2] += ene_half;
        } // END atom2 not excluded from atom1
      } // END atom1 or atom2 is selected
    } // END inner loop over atoms
  } // END outer loop over atoms
}
} // END namespace Energy
} // END namespace Cpptraj 
#endif

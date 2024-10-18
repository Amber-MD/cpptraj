#ifndef INC_ENERGY_ENE_NONBOND_H
#define INC_ENERGY_ENE_NONBOND_H
#include "Ene_LJ_6_12.h"
namespace Cpptraj {
namespace Energy {
/// Calculate simple nonbond energy with Coulomb and LJ 6-12
/** NOTE: To make this file more lightweight, there are no includes for
  *       Frame, AtomMask, etc. It is assumed that before this file is
  *       included there will be at least '#include "Topology.h" and
  *       'include <cmath>' (for the sqrt).
  */
template <typename T>
void Ene_Nonbond(Frame const& fIn, Topology const& tIn, AtomMask const& mask,
                 ExclusionArray const& Excluded, T const& QFAC, T& Eelec, T& Evdw)
{
  Evdw = 0;
  Eelec = 0;
  int idx1;
# ifdef _OPENMP
# pragma omp parallel private(idx1) reduction(+ : Eelec, Evdw)
  {
# pragma omp for
# endif
  for (idx1 = 0; idx1 < mask.Nselected(); idx1++)
  {
    int atom1 = mask[idx1];
    // Set up coord for this atom
    const double* crd1 = fIn.XYZ( atom1 );
    // Set up exclusion list for this atom
    // TODO refactor inner loop to be more like StructureCheck, more efficient.
    ExclusionArray::ExListType::const_iterator excluded_idx = Excluded[idx1].begin();
    for (int idx2 = idx1 + 1; idx2 < mask.Nselected(); idx2++)
    {
      int atom2 = mask[idx2];
      // Advance excluded list up to current selected atom
      while (excluded_idx != Excluded[idx1].end() && *excluded_idx < idx2) ++excluded_idx;
      // If atom is excluded, just increment to next excluded atom.
      if (excluded_idx != Excluded[idx1].end() && idx2 == *excluded_idx)
        ++excluded_idx;
      else {
        const double* crd2 = fIn.XYZ( atom2 );
        T dx = crd1[0] - crd2[0];
        T dy = crd1[1] - crd2[1];
        T dz = crd1[2] - crd2[2];
        T rij2 = (dx*dx) + (dy*dy) + (dz*dz);
        // VDW
        NonbondType const& LJ = tIn.GetLJparam(atom1, atom2);
        T e_vdw = Ene_LJ_6_12<T>( rij2, LJ.A(), LJ.B() );
        Evdw += e_vdw;
        // Coulomb
        T qiqj = QFAC * tIn[atom1].Charge() * tIn[atom2].Charge();
        T rij = sqrt( rij2 );
        T e_elec = qiqj / rij;
        Eelec += e_elec;
#       ifdef DEBUG_ENERGY
        mprintf("\tEVDW  %4i -- %4i: A=  %12.5e  B=  %12.5e  r2= %12.5f  E= %12.5e\n",
                atom1+1, atom2+1, LJ.A(), LJ.B(), rij2, e_vdw);
        mprintf("\tEELEC %4i -- %4i: q1= %12.5e  q2= %12.5e  r=  %12.5f  E= %12.5e\n",
                atom1+1, atom2+1, tIn[atom1].Charge(), tIn[atom2].Charge(),
                rij, e_elec);
#       endif
      }
    }
  }
# ifdef _OPENMP
  } // END omp parallel
# endif
}
} // END namespace Energy
} // END namespace Cpptraj
#endif

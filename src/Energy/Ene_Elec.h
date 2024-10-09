#ifndef INC_ENERGY_ENE_ELEC_H
#define INC_ENERGY_ENE_ELEC_H
namespace Cpptraj {
namespace Energy {
template <typename T>
T Ene_Elec(Frame const& fIn, Topology const& tIn, AtomMask const& mask, ExclusionArray const& Excluded)
{
  T Eelec = 0.0;
  int idx1;
# ifdef _OPENMP
# pragma omp parallel private(idx1) reduction(+ : Eelec)
  {
# pragma omp for
# endif
  for (idx1 = 0; idx1 < mask.Nselected(); idx1++)
  {
    int atom1 = mask[idx1];
    // Set up coord for this atom
    const double* xyz1 = fIn.XYZ( atom1 );
    // Set up exclusion list for this atom
    // TODO refactor inner loop to be like StructureCheck
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
        const double* xyz2 = fIn.XYZ( atom2 );
        T x = xyz1[0] - xyz2[0];
        T y = xyz1[1] - xyz2[1];
        T z = xyz1[2] - xyz2[2];
        T rij2 = (x*x + y*y + z*z);
        T rij = sqrt(rij2);
        // Coulomb
        T qiqj = Constants::COULOMBFACTOR * tIn[atom1].Charge() * tIn[atom2].Charge();
        T e_elec = qiqj / rij;
        Eelec += e_elec;
#       ifdef DEBUG_ENERGY
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
  return Eelec;
}
}
}
#endif

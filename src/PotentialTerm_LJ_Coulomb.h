#ifndef INC_POTENTIALTERM_LJ_COULOMB_H
#define INC_POTENTIALTERM_LJ_COULOMB_H
#include "PotentialTerm.h"
#include "ExclusionArray.h"
#include <vector>
// Forward declares
class NonbondParmType;
class Atom;
/// Potential term for simple nonbonded LJ + Coulomb with no cutoff
class PotentialTerm_LJ_Coulomb : public PotentialTerm {
  public:
    PotentialTerm_LJ_Coulomb();
    int InitTerm(MdOpts const&);
    int SetupTerm(Topology const&, Box const&, CharMask const&, EnergyArray&);
    void CalcForce(Frame&, CharMask const&) const;
  private:
    typedef std::vector<int> Iarray;
    Iarray selectedAtoms_;           ///< Selected atoms.
    Iarray nonSelectedAtoms_;        ///< Non-selected atoms.
    std::vector<Atom> const* atoms_; ///< Pointer to Atoms
    NonbondParmType const* nonbond_; ///< Pointer to nonbond parameters.
    double* E_vdw_;                  ///< Pointer to VDW term of energy array.
    double* E_elec_;                 ///< Pointer to Coulomb term of energy array.
    double QFAC_;                    ///< Coulomb calculation prefactor
    double cutoff2_;                 ///< Interaction distance cutoff squared (Ang^2)
    int nExclude_;                   ///< Interactions to exclude when setting up exclusion array
    ExclusionArray Excluded_;        ///< Exclusion array
};
#endif

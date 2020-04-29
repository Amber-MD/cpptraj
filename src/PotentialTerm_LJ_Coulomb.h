#ifndef INC_POTENTIALTERM_LJ_COULOMB_H
#define INC_POTENTIALTERM_LJ_COULOMB_H
#include "PotentialTerm.h"
#include <vector>
// Forward declares
class NonbondParmType;
/// Potential term for simple nonbonded LJ + Coulomb with no cutoff
class PotentialTerm_LJ_Coulomb : public PotentialTerm {
  public:
    PotentialTerm_LJ_Coulomb();
    int SetupTerm(Topology const&, CharMask const&, EnergyArray&);
    void CalcForce(Frame&, CharMask const&) const;
  private:
    typedef std::vector<int> Iarray;
    Iarray selectedAtoms_;           ///< Selected atoms.
    Iarray nonSelectedAtoms_;        ///< Non-selected atoms.
    Iarray typeIndices_;             ///< Type index for each atom
    NonbondParmType const* nonbond_; ///< Pointer to nonbond parameters.
    double* Evdw_;                   ///< Pointer to VDW term of energy array.
    double* Eelec_;                  ///< Pointer to Coulomb term of energy array.
};
#endif

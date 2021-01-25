#ifndef INC_POTENTIALTERM_DIHEDRAL_H
#define INC_POTENTIALTERM_DIHEDRAL_H
#include "PotentialTerm.h"
#include "ParameterTypes.h"
/// Torsion term, truncated Fourier series
class PotentialTerm_Dihedral : public PotentialTerm {
  public:
    PotentialTerm_Dihedral() : PotentialTerm(DIHEDRAL), dihParm_(0), Edih_(0) {}

    int SetupTerm(Topology const&, Box const&, CharMask const&, EnergyArray&);
    void CalcForce(Frame&, CharMask const&) const;
  private:
    void addDihedrals(DihedralArray const&, CharMask const&);

    DihedralArray activeDihs_;          ///< Array of dihedrals selected by mask during setup
    DihedralParmArray const* dihParm_;  ///< Pointer to array containing dihedral parameters
    double* Edih_;                      ///< Pointer to dihedral term of energy array.
};
#endif

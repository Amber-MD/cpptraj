#ifndef INC_POTENTIALTERM_DIHEDRAL_H
#define INC_POTENTIALTERM_DIHEDRAL_H
#include "PotentialTerm.h"
#include "ParameterTypes.h"
class Atom;
/// Torsion term, truncated Fourier series
class PotentialTerm_Dihedral : public PotentialTerm {
  public:
    PotentialTerm_Dihedral();

    int SetupTerm(Topology const&, Box const&, CharMask const&, EnergyArray&);
    void CalcForce(Frame&, CharMask const&) const;
  private:
    void addDihedrals(DihedralArray const&, CharMask const&);

    DihedralArray activeDihs_;          ///< Array of dihedrals selected by mask during setup
    DihedralParmArray const* dihParm_;  ///< Pointer to array containing dihedral parameters
    NonbondParmType const* nonbond_;    ///< Pointer to nonbond parameters
    std::vector<Atom> const* atoms_;    ///< Pointer to unmodified atoms array
    double* Edih_;                      ///< Pointer to dihedral term of energy array.
    double* Enb14_;                     ///< Pointer to vdw 14 term of energy array.
    double* Eq14_;                      ///< Pointer to elec 14 term of energy array
};
#endif

#ifndef INC_POTENTIALTERM_BOND_H
#define INC_POTENTIALTERM_BOND_H
#include "PotentialTerm.h"
#include "ParameterTypes.h"
/// Simple Hooke's law bond term
class PotentialTerm_Bond : public PotentialTerm {
  public:
    PotentialTerm_Bond() : PotentialTerm(BOND), bondParm_(0), Ebond_(0) {}

    int SetupTerm(Topology const&, Box const&, CharMask const&, EnergyArray&);
    void CalcForce(Frame&, CharMask const&) const;
  private:
    void addBonds(BondArray const&, CharMask const&);

    BondArray activeBonds_;         ///< Array of bonds selected by mask during setup
    BondParmArray const* bondParm_; ///< Pointer to array containing bond parameters
    double* Ebond_;                 ///< Pointer to bond term of energy array.
};
#endif

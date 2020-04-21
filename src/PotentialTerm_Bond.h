#ifndef INC_POTENTIALTERM_BOND_H
#define INC_POTENTIALTERM_BOND_H
#include "PotentialTerm.h"
#include "ParameterTypes.h"
/// Simple Hooke's law bond term
class PotentialTerm_Bond : public PotentialTerm {
  public:
    PotentialTerm_Bond() : PotentialTerm(BOND), bondParm_(0) {}

    int SetupTerm(Topology const&, CharMask const&);
    void CalcForce(Frame&, CharMask const&) const;
  private:
    BondArray activeBonds_;
    BondParmArray const* bondParm_;
};
#endif

#ifndef INC_POTENTIALTERM_OPENMM_H
#define INC_POTENTIALTERM_OPENMM_H
#include "PotentialTerm.h"
/// Use OpenMM potential
class PotentialTerm_OpenMM : public PotentialTerm {
  public:
    PotentialTerm_OpenMM() : PotentialTerm(OPENMM) {}

    int SetupTerm(Topology const&, CharMask const&, EnergyArray&);
    void CalcForce(Frame&, CharMask const&) const;
  private:
#   ifdef HAS_OPENMM
    OpenMM::System* system_;
    OpenMM::Context* context_;
#   endif
};
#endif

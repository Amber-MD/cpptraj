#ifndef INC_POTENTIALTERM_OPENMM_H
#define INC_POTENTIALTERM_OPENMM_H
#include "PotentialTerm.h"
// Forward declares
#ifdef HAS_OPENMM
namespace OpenMM {
  class System;
  class Context;
}
#endif
/// Use OpenMM potential
class PotentialTerm_OpenMM : public PotentialTerm {
  public:
    PotentialTerm_OpenMM();
    ~PotentialTerm_OpenMM();

    int SetupTerm(Topology const&, Box const&, CharMask const&, EnergyArray&);
    void CalcForce(Frame&, CharMask const&) const;
  private:
#   ifdef HAS_OPENMM
    OpenMM::System* system_;
    OpenMM::Context* context_;
#   endif
};
#endif

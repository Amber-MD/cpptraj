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

    int InitTerm(MdOpts const&);
    int SetupTerm(Topology const&, Box const&, CharMask const&, EnergyArray&);
    void CalcForce(Frame&, CharMask const&) const;
    int RemovedDegreesOfFreedom() const { return n_removed_dof_; }
  private:
#   ifdef HAS_OPENMM
    int OpenMM_setup(Topology const&, Box const&, CharMask const&, EnergyArray&);

    OpenMM::System* system_;
    OpenMM::Context* context_;
#   endif

    double* ene_;
    double scaleEE_; ///< Electrostatic 1-4 scaling factor
    double scaleNB_; ///< Lennard-Jones (VDW) 1-4 scaling factor
    double cut_;     ///< Nonbond cutoff, in nm
    int n_removed_dof_; ///< # of degrees of freedom removed by constraints
    bool shakeH_;    ///< If true, constain bonds to H
    bool shakeHeavy_;  ///< If true, constain bonds to heavy atoms
};
#endif

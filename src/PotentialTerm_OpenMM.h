#ifndef INC_POTENTIALTERM_OPENMM_H
#define INC_POTENTIALTERM_OPENMM_H
#include "PotentialTerm.h"
#include "ParameterTypes.h"
// Forward declares
#ifdef HAS_OPENMM
namespace OpenMM {
  class System;
  class Context;
  class HarmonicBondForce;
  class HarmonicAngleForce;
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
    void AddBonds(OpenMM::HarmonicBondForce*, std::vector< std::pair<int,int> >&,
                  BondArray const&, BondParmArray const&, std::vector<int> const&);
    void AddAngles(OpenMM::HarmonicAngleForce*, AngleArray const&, AngleParmArray const&, 
                   std::vector<int> const&);
    int OpenMM_setup(Topology const&, Box const&, CharMask const&, EnergyArray&);

    OpenMM::System* system_;
    OpenMM::Context* context_;
    double* ene_;
#   endif
};
#endif

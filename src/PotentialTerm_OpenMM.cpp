#include "PotentialTerm_OpenMM.h"
#include "CpptrajStdio.h"
#ifdef HAS_OPENMM
# include "OpenMM.h"
#endif

PotentialTerm_OpenMM::PotentialTerm_OpenMM() :
  PotentialTerm(OPENMM)
#ifdef HAS_OPENMM
  ,system_(0)
  ,context_(0)
#endif
{}

PotentialTerm_OpenMM::~PotentialTerm_OpenMM() {
# ifdef HAS_OPENMM
  if (system_ != 0) delete system_;
  if (context_ != 0) delete context_;
# endif
}

int PotentialTerm_OpenMM::SetupTerm(Topology const& topIn, CharMask const& maskIn,
                                    EnergyArray& earrayIn)
{
# ifdef HAS_OPENMM
  mprintf("OpenMM setup.\n");
  return 0;
# else
  mprinterr("Error: CPPTRAJ was compiled without OpenMM support.\n");
  return 1;
# endif
}

void PotentialTerm_OpenMM::CalcForce(Frame& frameIn, CharMask const& maskIn) const
{
  return;
}
  

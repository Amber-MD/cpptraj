#include "PotentialTerm_OpenMM.h"
#include "CpptrajStdio.h"
#ifdef HAS_OPENMM
# include "OpenMM.h"
# include "Box.h"
# include "Topology.h"
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

int PotentialTerm_OpenMM::SetupTerm(Topology const& topIn, Box const& boxIn,
                                    CharMask const& maskIn, EnergyArray& earrayIn)
{
# ifdef HAS_OPENMM
  mprintf("OpenMM setup.\n");
  system_ = new OpenMM::System();
  OpenMM::NonbondedForce* nonbond = new OpenMM::NonbondedForce();
  system_->addForce( nonbond );
  OpenMM::HarmonicBondForce* bondStretch = new OpenMM::HarmonicBondForce();
  system_->addForce( bondStretch );

  // Do periodic boundary conditions if necessary.
  if (boxIn.Type() != Box::NOBOX) {
    nonbond->setNonbondedMethod(OpenMM::NonbondedForce::CutoffPeriodic);
    nonbond->setCutoffDistance( 0.8 ); // TODO allow args
    Matrix_3x3 ucell, recip;
    boxIn.ToRecip(ucell, recip);
    system_->setDefaultPeriodicBoxVectors(
      OpenMM::Vec3( ucell[0], ucell[1], ucell[2] ),
      OpenMM::Vec3( ucell[3], ucell[4], ucell[5] ),
      OpenMM::Vec3( ucell[6], ucell[7], ucell[8] ) );
  }

  // Add atoms to the system.
  for (int idx = 0; idx != topIn.Natom(); idx++)
  {
    system_->addParticle( topIn[idx].Mass() );
    if (topIn.Nonbond().HasNonbond()) {
      nonbond->addParticle(
        topIn[idx].Charge(),
        topIn.GetVDWradius(idx) * OpenMM::NmPerAngstrom * OpenMM::SigmaPerVdwRadius,
        topIn.GetVDWdepth(idx) * OpenMM::KJPerKcal );
    }
  }
      
# else
  mprinterr("Error: CPPTRAJ was compiled without OpenMM support.\n");
  return 1;
# endif
}

void PotentialTerm_OpenMM::CalcForce(Frame& frameIn, CharMask const& maskIn) const
{
  return;
}
  

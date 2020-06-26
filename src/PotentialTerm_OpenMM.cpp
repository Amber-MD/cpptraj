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

#ifdef HAS_OPENMM
void PotentialTerm_OpenMM::AddBonds(OpenMM::HarmonicBondForce* bondStretch,
                                    BondArray const& bonds, BondParmArray const& BP)
{
  for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
    bondStretch->addBond( bnd->A1(), bnd->A2(),
                          BP[bnd->Idx()].Req() * OpenMM::NmPerAngstrom,
                          BP[bnd->Idx()].Rk() * 2 * OpenMM::KJPerKcal
                            * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm );
}


/** This performs the actual openMM setup. */
int PotentialTerm_OpenMM::OpenMM_setup(Topology const& topIn, Box const& boxIn,
                                       CharMask const& maskIn, EnergyArray& earrayIn)
{
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
  // Add bonds
  // Note factor of 2 for stiffness below because Amber specifies the constant
  // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants 
  // it as used in the force term kx, with energy kx^2/2.
  AddBonds(bondStretch, topIn.Bonds(), topIn.BondParm());
  AddBonds(bondStretch, topIn.BondsH(), topIn.BondParm());

  return 0;
}
#endif
 
/** Set up openmm terms. This is the wrapper for try/catch. */
int PotentialTerm_OpenMM::SetupTerm(Topology const& topIn, Box const& boxIn,
                                    CharMask const& maskIn, EnergyArray& earrayIn)
{
  int err = 1;
# ifdef HAS_OPENMM
  try {
    err = OpenMM_setup(topIn, boxIn, maskIn, earrayIn);
  }

  catch(const std::exception& e) {
    printf("EXCEPTION: %s\n", e.what());
    err = 1;
  }
# else
  mprinterr("Error: CPPTRAJ was compiled without OpenMM support.\n");
# endif
  return err;
}

void PotentialTerm_OpenMM::CalcForce(Frame& frameIn, CharMask const& maskIn) const
{
  return;
}
  

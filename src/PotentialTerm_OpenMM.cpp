#include "PotentialTerm_OpenMM.h"
#include "CpptrajStdio.h"
#ifdef HAS_OPENMM
# include "OpenMM.h"
# include "Box.h"
# include "Topology.h"
# include "CharMask.h"
# include "Constants.h"
# include "EnergyArray.h"
# include "PotentialTerm_InitOpts.h"
#endif

/// CONSTRUCTOR
PotentialTerm_OpenMM::PotentialTerm_OpenMM() :
  PotentialTerm(OPENMM)
#ifdef HAS_OPENMM
  ,system_(0)
  ,context_(0)
  ,scaleEE_(1.0/1.2) // Amber default
  ,scaleNB_(1.0/2.0) // Amber default
  ,cut_(0.8)         // in nm, Amber default
#endif
{}

/// DESTRUCTOR
PotentialTerm_OpenMM::~PotentialTerm_OpenMM() {
# ifdef HAS_OPENMM
  if (system_ != 0) delete system_;
  if (context_ != 0) delete context_;
# endif
}

#ifdef HAS_OPENMM
/** Add BondArray to openmm */
void PotentialTerm_OpenMM::AddBonds(OpenMM::HarmonicBondForce* bondStretch,
                                    OpenMM::System* system,
                                    std::vector< std::pair<int,int> >& bondPairs,
                                    BondArray const& bonds, BondParmArray const& BP,
                                    std::vector<int> const& oldToNew,
                                    bool useConstraints)
{
  for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
  {
    int a1 = oldToNew[bnd->A1()];
    int a2 = oldToNew[bnd->A2()];
    if (a1 != -1 && a2 != -1)
    {
      bondPairs.push_back(std::make_pair(a1, a2));
      if (useConstraints) {
        system->addConstraint( a1, a2, BP[bnd->Idx()].Req() * OpenMM::NmPerAngstrom );
      } else {
        bondStretch->addBond( a1, a2,
                              BP[bnd->Idx()].Req() * OpenMM::NmPerAngstrom,
                              BP[bnd->Idx()].Rk() * 2 * OpenMM::KJPerKcal
                                * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm );
      }
    }
  }
}

/** Add AngleArray to openmm. */
void PotentialTerm_OpenMM::AddAngles(OpenMM::HarmonicAngleForce* angleStretch,
                                     AngleArray const& angles, AngleParmArray const& AP,
                                     std::vector<int> const& oldToNew)
{
  for (AngleArray::const_iterator ang = angles.begin(); ang != angles.end(); ++ang)
  {
    int a1 = oldToNew[ang->A1()];
    int a2 = oldToNew[ang->A2()];
    int a3 = oldToNew[ang->A3()];
    if (a1 != -1 && a2 != -1 && a3 != -1)
    {
      // CPPTRAJ angles are already in radians, no need for OpenMM::RadiansPerDegree
      angleStretch->addAngle( a1, a2, a3,
                            AP[ang->Idx()].Teq(),
                            AP[ang->Idx()].Tk() * 2 * OpenMM::KJPerKcal);
    }
  }
}

/** Add DihedralArray to openmm. */
void PotentialTerm_OpenMM::AddDihedrals(OpenMM::PeriodicTorsionForce* ptorsion,
                                     DihedralArray const& dihedrals, DihedralParmArray const& DP,
                                     std::vector<int> const& oldToNew)
{
  for (DihedralArray::const_iterator dih = dihedrals.begin(); dih != dihedrals.end(); ++dih)
  {
    int a1 = oldToNew[dih->A1()];
    int a2 = oldToNew[dih->A2()];
    int a3 = oldToNew[dih->A3()];
    int a4 = oldToNew[dih->A4()];
    if (a1 != -1 && a2 != -1 && a3 != -1 && a4 != -1)
    {
      // CPPTRAJ dihedrals are already in radians, no need for OpenMM::RadiansPerDegree
      ptorsion->addTorsion( a1, a2, a3, a4,
                            DP[dih->Idx()].Pn(),
                            DP[dih->Idx()].Phase(),
                            DP[dih->Idx()].Pk() * OpenMM::KJPerKcal);
    }
  }
}

/** This performs the actual openMM setup. */
int PotentialTerm_OpenMM::OpenMM_setup(Topology const& topIn, Box const& boxIn,
                                       CharMask const& maskIn, EnergyArray& earrayIn)
{
  mprintf("\tSetting up OpenMM.\n");
  system_ = new OpenMM::System();
  OpenMM::NonbondedForce* nonbond = new OpenMM::NonbondedForce();
  system_->addForce( nonbond );
  OpenMM::HarmonicBondForce* bondStretch = new OpenMM::HarmonicBondForce();
  system_->addForce( bondStretch );
  OpenMM::HarmonicAngleForce* angleStretch = new OpenMM::HarmonicAngleForce();
  system_->addForce( angleStretch );
  OpenMM::PeriodicTorsionForce* ptorsion = new OpenMM::PeriodicTorsionForce();
  system_->addForce( ptorsion );

  // Do periodic boundary conditions if necessary.
  if (boxIn.Type() != Box::NOBOX) {
    nonbond->setNonbondedMethod(OpenMM::NonbondedForce::CutoffPeriodic);
    nonbond->setCutoffDistance( cut_ ); // TODO allow args
    Matrix_3x3 ucell, recip;
    boxIn.ToRecip(ucell, recip);
    system_->setDefaultPeriodicBoxVectors(
      OpenMM::Vec3( ucell[0], ucell[1], ucell[2] ),
      OpenMM::Vec3( ucell[3], ucell[4], ucell[5] ),
      OpenMM::Vec3( ucell[6], ucell[7], ucell[8] ) );
  }

  // Add atoms to the system.
  std::vector<int> oldToNew(topIn.Natom(), -1);
  int newIdx = 0;
  for (int idx = 0; idx != topIn.Natom(); idx++)
  {
    if (maskIn.AtomInCharMask(idx)) {
      oldToNew[idx] = newIdx++;
      system_->addParticle( topIn[idx].Mass() );
      if (topIn.Nonbond().HasNonbond()) {
        nonbond->addParticle(
          topIn[idx].Charge(),
          topIn.GetVDWradius(idx) * OpenMM::NmPerAngstrom * OpenMM::SigmaPerVdwRadius,
          topIn.GetVDWdepth(idx) * OpenMM::KJPerKcal );
      }
    }
  }
  // Add bonds
  // Note factor of 2 for stiffness below because Amber specifies the constant
  // as it is used in the harmonic energy term kx^2 with force 2kx; OpenMM wants 
  // it as used in the force term kx, with energy kx^2/2.
  std::vector< std::pair<int,int> >   bondPairs;
  AddBonds(bondStretch, system_, bondPairs, topIn.Bonds(), topIn.BondParm(),  oldToNew, false);
  AddBonds(bondStretch, system_, bondPairs, topIn.BondsH(), topIn.BondParm(), oldToNew, false);

  // Add angles
  AddAngles(angleStretch, topIn.Angles(), topIn.AngleParm(), oldToNew);
  AddAngles(angleStretch, topIn.AnglesH(), topIn.AngleParm(), oldToNew);

  // Add dihedrals
  AddDihedrals(ptorsion, topIn.Dihedrals(), topIn.DihedralParm(), oldToNew);
  AddDihedrals(ptorsion, topIn.DihedralsH(), topIn.DihedralParm(), oldToNew);

  // Populate nonbonded exclusions TODO make args
  nonbond->createExceptionsFromBonds(bondPairs, scaleEE_, scaleNB_);

  // Set up integrator and context.
  OpenMM::Integrator* integrator = new OpenMM::VerletIntegrator( 0.001 ); // TODO allow ars
  context_ = new OpenMM::Context( *system_, *integrator );

  std::string platformName = context_->getPlatform().getName();
  mprintf("\tOpenMM Platform: %s\n", platformName.c_str());

  // Set up energy term
  ene_ = earrayIn.AddType( EnergyArray::E_OPENMM );

  return 0;
}
#endif

/** Init openmm options. */
int PotentialTerm_OpenMM::InitTerm(InitOpts const& opts) {
# ifdef HAS_OPENMM
  scaleEE_ = opts.ScaleEE();
  scaleNB_ = opts.ScaleNB();
  if (opts.CutEE() != opts.CutNB()) {
    mprinterr("Error: Elec. cut %g != VDW cut %g; not yet supported.\n");
    return 1;
  }
  cut_ = opts.CutEE();
# else
  mprinterr("Error: CPPTRAJ was compiled without OpenMM support.\n");
  return 1;
# endif
}

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

/** Calculate force from OpenMM */
void PotentialTerm_OpenMM::CalcForce(Frame& frameIn, CharMask const& maskIn) const
{
# ifdef HAS_OPENMM
  // Set positions, convert from Ang to nm
  std::vector<OpenMM::Vec3> posInNm;
  posInNm.reserve(maskIn.Nselected());
  for (int at = 0; at != frameIn.Natom(); at++) {
    if (maskIn.AtomInCharMask(at)) {
      const double* xyz = frameIn.XYZ(at);
      posInNm.push_back( OpenMM::Vec3(xyz[0]*OpenMM::NmPerAngstrom,
                                      xyz[1]*OpenMM::NmPerAngstrom,
                                      xyz[2]*OpenMM::NmPerAngstrom) );
    }
  }
  context_->setPositions(posInNm);
  // Do a single minimization step
  //OpenMM::LocalEnergyMinimizer min;
  //min.minimize(*context_, 10.0, 1);
  // Get the results
  const OpenMM::State state = context_->getState(OpenMM::State::Forces |
                                                 OpenMM::State::Energy, true);
  //  timeInPs = state.getTime(); // OpenMM time is in ps already
  // Convert to kcal/mol from kJ/mol
  *ene_ = state.getPotentialEnergy() * Constants::J_TO_CAL;

  // Copy OpenMM forces into output array and change units from nm/kJ to Angstroms/kcal.
  const std::vector<OpenMM::Vec3>& ommForces = state.getForces();
  double* fptr = frameIn.fAddress();

  int oi = 0; // Index into ommForces
  for (int at = 0; at != frameIn.Natom(); at++, fptr += 3)
  {
    if (maskIn.AtomInCharMask(at)) {
      for (int j = 0; j < 3; j++)
        fptr[j] += ommForces[oi][j] * Constants::GMX_FRC_TO_AMBER;
       oi++;
     }
  }

# endif
  return;
}

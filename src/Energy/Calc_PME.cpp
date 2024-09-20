#include "Calc_PME.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"
#include "../Topology.h"

using namespace Cpptraj::Energy;

Calc_PME::Calc_PME() :
  Recip_(PME_Recip::COULOMB)
{}

/** Set up PME parameters. */
int Calc_PME::Init(Box const& boxIn, EwaldOptions const& pmeOpts, int debugIn)
{
  if (NBengine_.ModifyEwaldParams().InitEwald(boxIn, pmeOpts, debugIn)) {
    mprinterr("Error: PME calculation init failed.\n");
    return 1;
  }
  if (pairList_.InitPairList(NBengine_.EwaldParams().Cutoff(), pmeOpts.SkinNB(), debugIn))
    return 1;
  if (pairList_.SetupPairList( boxIn ))
    return 1;
  VDW_LR_.SetDebug( debugIn );
  Recip_.SetDebug( debugIn );

  return 0;
}

/** Setup PME calculation. */
int Calc_PME::Setup(Topology const& topIn, AtomMask const& maskIn) {
  if (NBengine_.ModifyEwaldParams().SetupEwald(topIn, maskIn)) {
    mprinterr("Error: PME calculation setup failed.\n");
    return 1;
  }
  if (VDW_LR_.Setup_VDW_Correction( topIn, maskIn )) {
    mprinterr("Error: PME calculation long range VDW correction setup failed.\n");
    return 1;
  }
  // Setup exclusion list
  // Use distance of 4 (up to dihedrals)
  if (Excluded_.SetupExcluded(topIn.Atoms(), maskIn, 4,
                              ExclusionArray::EXCLUDE_SELF,
                              ExclusionArray::FULL))
  {
    mprinterr("Error: Could not set up exclusion list for PME calculation.\n");
    return 1;
  }

  return 0;
}

/** Calculate full nonbonded energy with PME */
int Calc_PME::CalcNonbondEnergy(Frame const& frameIn, AtomMask const& maskIn,
                                double& e_elec, double& e_vdw)
{
  t_total_.Start();
  double volume = frameIn.BoxCrd().CellVolume();
  double e_self = Self( volume );
  double e_vdw_lr_correction;

  int retVal = pairList_.CreatePairList(frameIn, frameIn.BoxCrd().UnitCell(), frameIn.BoxCrd().FracCell(), maskIn);
  if (retVal != 0) {
    mprinterr("Error: Grid setup failed.\n");
    return 1;
  }

  // TODO make more efficient
  int idx = 0;
  coordsD_.clear();
  for (AtomMask::const_iterator atm = maskIn.begin(); atm != maskIn.end(); ++atm, ++idx) {
    const double* XYZ = frameIn.XYZ( *atm );
    coordsD_.push_back( XYZ[0] );
    coordsD_.push_back( XYZ[1] );
    coordsD_.push_back( XYZ[2] );
  }

//  MapCoords(frameIn, ucell, recip, maskIn);
  double e_recip = Recip_ParticleMesh( frameIn.BoxCrd() );

  // TODO branch
  double e_vdw6self, e_vdw6recip;
  if (lw_coeff_ > 0.0) {
    e_vdw6self = Self6();
    e_vdw6recip = LJ_Recip_ParticleMesh( frameIn.BoxCrd() );
    if (debug_ > 0) {
      mprintf("DEBUG: e_vdw6self = %16.8f\n", e_vdw6self);
      mprintf("DEBUG: Evdwrecip = %16.8f\n", e_vdw6recip);
    }
    e_vdw_lr_correction = 0.0;
  } else {
    e_vdw6self = 0.0;
    e_vdw6recip = 0.0;
    e_vdw_lr_correction = Vdw_Correction( volume );
  }

  e_vdw = 0.0;
  double e_adjust = 0.0;
  double e_direct = Direct( pairList_, e_vdw, e_adjust );
  if (debug_ > 0) {
    mprintf("DEBUG: Nonbond energy components:\n");
    mprintf("     Evdw                   = %24.12f\n", e_vdw + e_vdw_lr_correction + e_vdw6self + e_vdw6recip);
    mprintf("     Ecoulomb               = %24.12f\n", e_self + e_recip + e_direct + e_adjust);
    mprintf("\n");
    mprintf("     E electrostatic (self) = %24.12f\n", e_self);
    mprintf("                     (rec)  = %24.12f\n", e_recip);
    mprintf("                     (dir)  = %24.12f\n", e_direct);
    mprintf("                     (adj)  = %24.12f\n", e_adjust);
    mprintf("     E vanDerWaals   (dir)  = %24.12f\n", e_vdw);
    mprintf("                     (LR)   = %24.12f\n", e_vdw_lr_correction);
    mprintf("                     (6slf) = %24.12f\n", e_vdw6self);
    mprintf("                     (6rcp) = %24.12f\n", e_vdw6recip);
  }
  e_vdw += (e_vdw_lr_correction + e_vdw6self + e_vdw6recip);
  t_total_.Stop();
  e_elec = e_self + e_recip + e_direct + e_adjust;
  return 0;
}



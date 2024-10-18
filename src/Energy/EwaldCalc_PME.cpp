#include "EwaldCalc_PME.h"
#ifdef LIBPME
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"
#include "../Frame.h"
#include "../PairListTemplate.h"

using namespace Cpptraj::Energy;

EwaldCalc_PME::EwaldCalc_PME() :
  Recip_(PME_Recip::COULOMB)
{}

/** Set up PME parameters. */
int EwaldCalc_PME::Init(Box const& boxIn, EwaldOptions const& pmeOpts, int debugIn)
{
  // Sanity check
  if (!pmeOpts.IsPmeType()) {
    mprinterr("Internal Error: EwaldCalc_PME::Init(): Options were not set up for PME.\n");
    return 1;
  }
  mprintf("\tParticle Mesh Ewald params:\n");
  if (NBengine_.ModifyEwaldParams().InitEwald(boxIn, pmeOpts, debugIn)) {
    mprinterr("Error: PME calculation init failed.\n");
    return 1;
  }
  VDW_LR_.SetDebug( debugIn );
  if (Recip_.InitRecip(pmeOpts, debugIn)) {
    mprinterr("Error: PME recip init failed.\n");
    return 1;
  }

  return 0;
}

/** Setup PME calculation. */
int EwaldCalc_PME::Setup(Topology const& topIn, AtomMask const& maskIn) {
  if (NBengine_.ModifyEwaldParams().SetupEwald(topIn, maskIn)) {
    mprinterr("Error: PME calculation setup failed.\n");
    return 1;
  }
  if (VDW_LR_.Setup_VDW_Correction( topIn, maskIn )) {
    mprinterr("Error: PME calculation long range VDW correction setup failed.\n");
    return 1;
  }

  return 0;
}

/** Calculate full nonbonded energy with PME */
int EwaldCalc_PME::CalcNonbondEnergy(Frame const& frameIn, AtomMask const& maskIn,
                                     PairList const& pairList_, ExclusionArray const& Excluded_,
                                     double& e_elec, double& e_vdw)
{
  t_total_.Start();
  double volume = frameIn.BoxCrd().CellVolume();
  double e_self = NBengine_.EwaldParams().SelfEnergy( volume );

  // TODO make more efficient
  NBengine_.ModifyEwaldParams().FillRecipCoords( frameIn, maskIn );

  // FIXME helPME requires coords and charge arrays to be non-const
  double e_recip = Recip_.Recip_ParticleMesh( NBengine_.ModifyEwaldParams().SelectedCoords(),
                                              frameIn.BoxCrd(),
                                              NBengine_.ModifyEwaldParams().SelectedCharges(),
                                              NBengine_.EwaldParams().EwaldCoeff()
                                            );

  double e_vdw_lr_correction = VDW_LR_.Vdw_Correction( NBengine_.EwaldParams().Cutoff(), volume );

  t_direct_.Start();
  Cpptraj::PairListTemplate<double>(pairList_, Excluded_,
                                    NBengine_.EwaldParams().Cut2(), NBengine_);
  t_direct_.Stop();

  if (NBengine_.EwaldParams().Debug() > 0) {
    mprintf("DEBUG: Nonbond energy components:\n");
    mprintf("     Evdw                   = %24.12f\n", NBengine_.Evdw() + e_vdw_lr_correction );
    mprintf("     Ecoulomb               = %24.12f\n", e_self + e_recip +
                                                       NBengine_.Eelec() +
                                                       NBengine_.Eadjust());
    mprintf("\n");
    mprintf("     E electrostatic (self) = %24.12f\n", e_self);
    mprintf("                     (rec)  = %24.12f\n", e_recip);
    mprintf("                     (dir)  = %24.12f\n", NBengine_.Eelec());
    mprintf("                     (adj)  = %24.12f\n", NBengine_.Eadjust());
    mprintf("     E vanDerWaals   (dir)  = %24.12f\n", NBengine_.Evdw());
    mprintf("                     (LR)   = %24.12f\n", e_vdw_lr_correction);
  }
  e_vdw = NBengine_.Evdw() + e_vdw_lr_correction;
  e_elec = e_self + e_recip + NBengine_.Eelec() + NBengine_.Eadjust();
  t_total_.Stop();
  return 0;
}

void EwaldCalc_PME::Timing(double total) const {
  t_total_.WriteTiming(1,  "  PME Total:", total);
  Recip_.Timing_Total().WriteTiming(2,  "Recip:     ", t_total_.Total());
  t_direct_.WriteTiming(2, "Direct:    ", t_total_.Total());
}
#endif /* LIBPME */

#include "EwaldCalc_Regular.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"
#include "../Frame.h"
#include "../PairListTemplate.h"

using namespace Cpptraj::Energy;

EwaldCalc_Regular::EwaldCalc_Regular()
{}

/** Set up regular Ewald parameters. */
int EwaldCalc_Regular::Init(Box const& boxIn, EwaldOptions const& pmeOpts, int debugIn)
{
  // Sanity check
  if (pmeOpts.IsPmeType()) {
    mprinterr("Internal Error: EwaldCalc_Regular::Init(): Options were set up for PME only.\n");
    return 1;
  }
  mprintf("\tRegular Ewald params:\n");
  if (NBengine_.ModifyEwaldParams().InitEwald(boxIn, pmeOpts, debugIn)) {
    mprinterr("Error: Ewald calculation init failed.\n");
    return 1;
  }
  VDW_LR_.SetDebug( debugIn );
  if (Recip_.InitRecip(pmeOpts, NBengine_.EwaldParams().EwaldCoeff(), boxIn, debugIn)) {
    mprinterr("Error: Ewald recip init failed.\n");
    return 1;
  }

  return 0;
}

/** Setup PME calculation. */
int EwaldCalc_Regular::Setup(Topology const& topIn, AtomMask const& maskIn) {
  if (NBengine_.ModifyEwaldParams().SetupEwald(topIn, maskIn)) {
    mprinterr("Error: Ewald calculation setup failed.\n");
    return 1;
  }
  if (VDW_LR_.Setup_VDW_Correction( topIn, maskIn )) {
    mprinterr("Error: Ewald calculation long range VDW correction setup failed.\n");
    return 1;
  }
  if (Recip_.SetupRecip( maskIn.Nselected() )) {
    mprinterr("Error: Ewald calculation recip setup failed.\n");
    return 1;
  }

  return 0;
}

/** Calculate full nonbonded energy with regular Ewald */
int EwaldCalc_Regular::CalcNonbondEnergy(Frame const& frameIn, AtomMask const& maskIn,
                                     PairList const& pairList_, ExclusionArray const& Excluded_,
                                     double& e_elec, double& e_vdw)
{
  t_total_.Start();
  double volume = frameIn.BoxCrd().CellVolume();
  double e_self = NBengine_.EwaldParams().SelfEnergy( volume );

  // TODO make more efficient
  //NBengine_.ModifyEwaldParams().FillRecipCoords( frameIn, maskIn );

  double e_recip = Recip_.Recip_Regular( frameIn.BoxCrd().FracCell(),
                                         volume,
                                         pairList_.FracCoords(),
                                         NBengine_.EwaldParams().Charge() );

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

void EwaldCalc_Regular::Timing(double total) const {
  t_total_.WriteTiming(1,  "  Ewald Total:", total);
  Recip_.PrintTiming(t_total_.Total());
  t_direct_.WriteTiming(2, "Direct:    ", t_total_.Total());
}

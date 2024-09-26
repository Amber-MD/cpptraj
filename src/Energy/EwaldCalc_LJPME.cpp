#include "EwaldCalc_LJPME.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"
#include "../PairListTemplate.h"
#include "../Topology.h"

using namespace Cpptraj::Energy;

EwaldCalc_LJPME::EwaldCalc_LJPME() :
  Recip_(PME_Recip::COULOMB),
  LJrecip_(PME_Recip::LJ)
{}

/** Set up LJPME parameters. */
int EwaldCalc_LJPME::Init(Box const& boxIn, EwaldOptions const& pmeOpts, int debugIn)
{
  if (NBengine_.ModifyEwaldParams().InitEwald(boxIn, pmeOpts, debugIn)) {
    mprinterr("Error: LJPME calculation init failed.\n");
    return 1;
  }
  if (pairList_.InitPairList(NBengine_.EwaldParams().Cutoff(), pmeOpts.SkinNB(), debugIn))
    return 1;
  if (pairList_.SetupPairList( boxIn ))
    return 1;
  Recip_.SetDebug( debugIn );
  LJrecip_.SetDebug( debugIn );

  return 0;
}

/** Setup LJPME calculation. */
int EwaldCalc_LJPME::Setup(Topology const& topIn, AtomMask const& maskIn) {
  if (NBengine_.ModifyEwaldParams().SetupEwald(topIn, maskIn)) {
    mprinterr("Error: LJPME calculation setup failed.\n");
    return 1;
  }
  // Setup exclusion list
  // Use distance of 4 (up to dihedrals)
  if (Excluded_.SetupExcluded(topIn.Atoms(), maskIn, 4,
                              ExclusionArray::EXCLUDE_SELF,
                              ExclusionArray::FULL))
  {
    mprinterr("Error: Could not set up exclusion list for LJPME calculation.\n");
    return 1;
  }

  return 0;
}

/** Calculate full nonbonded energy with LJPME */
int EwaldCalc_LJPME::CalcNonbondEnergy(Frame const& frameIn, AtomMask const& maskIn,
                                double& e_elec, double& e_vdw)
{
  t_total_.Start();
  double volume = frameIn.BoxCrd().CellVolume();
  double e_self = NBengine_.EwaldParams().SelfEnergy( volume );

  int retVal = pairList_.CreatePairList(frameIn, frameIn.BoxCrd().UnitCell(),
                                        frameIn.BoxCrd().FracCell(), maskIn);
  if (retVal != 0) {
    mprinterr("Error: Pairlist creation failed for LJPME calc.\n");
    return 1;
  }

  // TODO make more efficient
  NBengine_.ModifyEwaldParams().FillRecipCoords( frameIn, maskIn );

  //  MapCoords(frameIn, ucell, recip, maskIn);
  // FIXME helPME requires coords and charge arrays to be non-const
  double e_recip = Recip_.Recip_ParticleMesh( NBengine_.ModifyEwaldParams().SelectedCoords(),
                                              frameIn.BoxCrd(),
                                              NBengine_.ModifyEwaldParams().SelectedCharges(),
                                              NBengine_.EwaldParams().NFFT(),
                                              NBengine_.EwaldParams().EwaldCoeff(),
                                              NBengine_.EwaldParams().Order()
                                            );
  double e_vdw6recip = LJrecip_.Recip_ParticleMesh( NBengine_.ModifyEwaldParams().SelectedCoords(),
                                                    frameIn.BoxCrd(),
                                                    NBengine_.ModifyEwaldParams().SelectedC6params(),
                                                    NBengine_.EwaldParams().NFFT(),
                                                    NBengine_.EwaldParams().LJ_EwaldCoeff(),
                                                    NBengine_.EwaldParams().Order()
                                                  );
  if (NBengine_.EwaldParams().Debug() > 0) {
    mprintf("DEBUG: e_vdw6self = %16.8f\n", NBengine_.EwaldParams().Self6());
    mprintf("DEBUG: Evdwrecip = %16.8f\n", e_vdw6recip);
  }

  t_direct_.Start();
  Cpptraj::PairListTemplate<double>(pairList_, Excluded_,
                                    NBengine_.EwaldParams().Cut2(), NBengine_);
  t_direct_.Stop();

  if (NBengine_.EwaldParams().Debug() > 0) {
    mprintf("DEBUG: Nonbond energy components:\n");
    mprintf("     Evdw                   = %24.12f\n", NBengine_.Evdw() +
                                                       NBengine_.EwaldParams().Self6() +
                                                       e_vdw6recip);
    mprintf("     Ecoulomb               = %24.12f\n", e_self + e_recip +
                                                       NBengine_.Eelec() +
                                                       NBengine_.Eadjust());
    mprintf("\n");
    mprintf("     E electrostatic (self) = %24.12f\n", e_self);
    mprintf("                     (rec)  = %24.12f\n", e_recip);
    mprintf("                     (dir)  = %24.12f\n", NBengine_.Eelec());
    mprintf("                     (adj)  = %24.12f\n", NBengine_.Eadjust());
    mprintf("     E vanDerWaals   (dir)  = %24.12f\n", NBengine_.Evdw());
    mprintf("                     (6slf) = %24.12f\n", NBengine_.EwaldParams().Self6());
    mprintf("                     (6rcp) = %24.12f\n", e_vdw6recip);
  }
  e_vdw = NBengine_.Evdw() + NBengine_.EwaldParams().Self6() + e_vdw6recip;
  e_elec = e_self + e_recip + NBengine_.Eelec() + NBengine_.Eadjust();
  t_total_.Stop();
  return 0;
}

void EwaldCalc_LJPME::Timing(double total) const {
  t_total_.WriteTiming(1,  "  LJPME Total:", total);
  Recip_.Timing_Total().WriteTiming(2,  "Recip:     ", t_total_.Total());
  //Recip_.Timing_Calc().WriteTiming(3,  "Recip. Calc   :", Recip_.Timing_Total().Total());
  LJrecip_.Timing_Total().WriteTiming(2,"LJRecip:   ", t_total_.Total());
  //LJrecip_.Timing_Calc().WriteTiming(3,"LJ Recip. Calc:", LJrecip_.Timing_Total().Total());
  t_direct_.WriteTiming(2, "Direct:    ", t_total_.Total());

  pairList_.Timing(total);
}
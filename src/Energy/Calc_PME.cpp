#include "Calc_PME.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"
#include "../PairListTemplate.h"
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
  double e_self = NBengine_.EwaldParams().SelfEnergy( volume );

  int retVal = pairList_.CreatePairList(frameIn, frameIn.BoxCrd().UnitCell(),
                                        frameIn.BoxCrd().FracCell(), maskIn);
  if (retVal != 0) {
    mprinterr("Error: Pairlist creation failed for PME calc.\n");
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

  // TODO branch
  //double e_vdw6self, e_vdw6recip;
  //if (lw_coeff_ > 0.0) {
  //  e_vdw6self = Self6();
  //  e_vdw6recip = LJ_Recip_ParticleMesh( frameIn.BoxCrd() );
  //  if (debug_ > 0) {
  //    mprintf("DEBUG: e_vdw6self = %16.8f\n", e_vdw6self);
  //    mprintf("DEBUG: Evdwrecip = %16.8f\n", e_vdw6recip);
  //  }
  //  e_vdw_lr_correction = 0.0;
  //} else {
  //  e_vdw6self = 0.0;
  //  e_vdw6recip = 0.0;
    double e_vdw_lr_correction = VDW_LR_.Vdw_Correction( NBengine_.EwaldParams().Cutoff(), volume );
  //}
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

void Calc_PME::Timing(double total) const {
  t_total_.WriteTiming(1,  "  PME Total:", total);
  //t_self_.WriteTiming(2,   "Self:      ", t_total_.Total());
  Recip_.Timing().WriteTiming(2,  "Recip:     ", t_total_.Total());
//  if (t_trig_tables_.Total() > 0.0)
//    t_trig_tables_.WriteTiming(3, "Calc trig tables:", t_recip_.Total());
  t_direct_.WriteTiming(2, "Direct:    ", t_total_.Total());
//# ifndef _OPENMP
//  t_erfc_.WriteTiming(3,  "ERFC:  ", t_direct_.Total());
//  t_adjust_.WriteTiming(3,"Adjust:", t_direct_.Total());
//# endif
  pairList_.Timing(total);
}

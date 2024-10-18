#include "EwaldCalc_Decomp_PME.h"
#ifdef LIBPME
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"
#include "../Frame.h"
#include "../PairListTemplate.h"

using namespace Cpptraj::Energy;

EwaldCalc_Decomp_PME::EwaldCalc_Decomp_PME() :
  Recip_(PME_Recip::COULOMB)
{}

/** Set up PME parameters. */
int EwaldCalc_Decomp_PME::Init(Box const& boxIn, EwaldOptions const& pmeOpts, int debugIn)
{
  // Sanity check
  if (!pmeOpts.IsPmeType()) {
    mprinterr("Internal Error: EwaldCalc_Decomp_PME::Init(): Options were not set up for PME.\n");
    return 1;
  }
  mprintf("\tDecomposed Particle Mesh Ewald params:\n");
  if (NBengine_.ModifyEwaldParams().InitEwald(boxIn, pmeOpts, debugIn)) {
    mprinterr("Error: Decomposable PME calculation init failed.\n");
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
int EwaldCalc_Decomp_PME::Setup(Topology const& topIn, AtomMask const& maskIn) {
  if (NBengine_.ModifyEwaldParams().SetupEwald(topIn, maskIn)) {
    mprinterr("Error: PME calculation setup failed.\n");
    return 1;
  }
  if (VDW_LR_.Setup_VDW_Correction( topIn, maskIn )) {
    mprinterr("Error: PME calculation long range VDW correction setup failed.\n");
    return 1;
  }
  // TODO reserve atom_elec and atom_vdw?

  return 0;
}

/** Calculate full decompoesd nonbonded energy with PME */
int EwaldCalc_Decomp_PME::CalcNonbondEnergy(Frame const& frameIn, AtomMask const& maskIn,
                                            PairList const& pairList_, ExclusionArray const& Excluded_,
                                            double& e_elec, double& e_vdw)
{
  t_total_.Start();
  double volume = frameIn.BoxCrd().CellVolume();
  Darray atom_self;
  double e_self = NBengine_.EwaldParams().DecomposedSelfEnergy( atom_self, volume );
# ifdef CPPTRAJ_DEBUG_ENEDECOMP
  mprintf("DEBUG: Total self energy: %f\n", e_self);
  mprintf("DEBUG: Sum of self array: %f\n", sumArray(atom_self));
# endif
  // TODO make more efficient
  NBengine_.ModifyEwaldParams().FillRecipCoords( frameIn, maskIn );

  Darray atom_recip;
  // FIXME helPME requires coords and charge arrays to be non-const
  double e_recip = Recip_.Recip_Decomp( atom_recip,
                                        NBengine_.ModifyEwaldParams().SelectedCoords(),
                                        frameIn.BoxCrd(),
                                        NBengine_.ModifyEwaldParams().SelectedCharges(),
                                        NBengine_.EwaldParams().EwaldCoeff()
                                      );
# ifdef CPPTRAJ_DEBUG_ENEDECOMP
  mprintf("DEBUG: Recip energy      : %f\n", e_recip);
  mprintf("DEBUG: Sum of recip array: %f\n", sumArray(atom_recip));
# endif
  // TODO distribute
  Darray atom_vdwlr;
  double e_vdw_lr_correction = VDW_LR_.Vdw_Decomp_Correction( atom_vdwlr,
                                                       NBengine_.EwaldParams().Cutoff(), volume );
# ifdef CPPTRAJ_DEBUG_ENEDECOMP
  mprintf("DEBUG: VDW correction       : %f\n", e_vdw_lr_correction);
  mprintf("DEBUG: Sum of VDW correction: %f\n", sumArray(atom_vdwlr));
# endif
  t_direct_.Start();
  Cpptraj::PairListTemplate<double>(pairList_, Excluded_,
                                    NBengine_.EwaldParams().Cut2(), NBengine_);
  t_direct_.Stop();
# ifdef CPPTRAJ_DEBUG_ENEDECOMP
  mprintf("DEBUG: Direct Elec. energy : %f\n", NBengine_.Eelec());
  mprintf("DEBUG: Sum of elec. energy : %f\n", sumArray(NBengine_.Eatom_Elec()));
  mprintf("DEBUG: Direct VDW energy   : %f\n", NBengine_.Evdw());
  mprintf("DEBUG: Sum of VDW energy   : %f\n", sumArray(NBengine_.Eatom_EVDW()));
  mprintf("DEBUG: Direct Adjust energy: %f\n", NBengine_.Eadjust());
  mprintf("DEBUG: Sum of Adjust energy: %f\n", sumArray(NBengine_.Eatom_EAdjust()));
# endif
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
  // TODO preallocate?
  atom_elec_.resize( NBengine_.Eatom_Elec().size() );
  atom_vdw_.resize(  NBengine_.Eatom_EVDW().size() );
  for (unsigned int idx = 0; idx != atom_elec_.size(); idx++)
  {
    atom_elec_[idx] = atom_self[idx] + atom_recip[idx] + NBengine_.Eatom_Elec()[idx] + NBengine_.Eatom_EAdjust()[idx];
    atom_vdw_[idx]  = NBengine_.Eatom_EVDW()[idx] + atom_vdwlr[idx];
  }
  t_total_.Stop();
  return 0;
}

void EwaldCalc_Decomp_PME::Timing(double total) const {
  t_total_.WriteTiming(1,  "  PME decomp Total:", total);
  Recip_.Timing_Total().WriteTiming(2,  "Recip:     ", t_total_.Total());
  t_direct_.WriteTiming(2, "Direct:    ", t_total_.Total());
}
#endif /* LIBPME */

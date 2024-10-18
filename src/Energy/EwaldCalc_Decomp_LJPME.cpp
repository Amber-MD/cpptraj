#include "EwaldCalc_Decomp_LJPME.h"
#ifdef LIBPME
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"
#include "../Frame.h"
#include "../PairListTemplate.h"

using namespace Cpptraj::Energy;

EwaldCalc_Decomp_LJPME::EwaldCalc_Decomp_LJPME() :
  Recip_(PME_Recip::COULOMB),
  LJrecip_(PME_Recip::LJ)
{}

/** Set up LJPME parameters. */
int EwaldCalc_Decomp_LJPME::Init(Box const& boxIn, EwaldOptions const& pmeOpts, int debugIn)
{
  // Sanity check
  if (!pmeOpts.IsPmeType()) {
    mprinterr("Internal Error: EwaldCalc_Decomp_LJPME::Init(): Options were not set up for PME.\n");
    return 1;
  }
  mprintf("\tDecomposable Particle Mesh Ewald (LJPME) params:\n");
  if (NBengine_.ModifyEwaldParams().InitEwald(boxIn, pmeOpts, debugIn)) {
    mprinterr("Error: Decomposable LJPME calculation init failed.\n");
    return 1;
  }
  if (Recip_.InitRecip(pmeOpts, debugIn)) {
    mprinterr("Error: Decomposable LJPME recip init failed.\n");
    return 1;
  }
  if (LJrecip_.InitRecip(pmeOpts, debugIn)) {
    mprinterr("Error: Decomposable LJPME LJ recip init failed.\n");
    return 1;
  }

  return 0;
}

/** Setup LJPME calculation. */
int EwaldCalc_Decomp_LJPME::Setup(Topology const& topIn, AtomMask const& maskIn) {
  if (NBengine_.ModifyEwaldParams().SetupEwald(topIn, maskIn)) {
    mprinterr("Error: LJPME calculation setup failed.\n");
    return 1;
  }
  NBengine_.ModifyEwaldParams().CalcDecomposedSelf6Energy();
  // TODO reserve atom_elec and atom_vdw?

  return 0;
}

/** Calculate full nonbonded energy with LJPME */
int EwaldCalc_Decomp_LJPME::CalcNonbondEnergy(Frame const& frameIn, AtomMask const& maskIn,
                                              PairList const& pairList_, ExclusionArray const& Excluded_,
                                              double& e_elec, double& e_vdw)
{
  t_total_.Start();
  double volume = frameIn.BoxCrd().CellVolume();
  // Do decomposed self
  Darray atom_self;
  double e_self = NBengine_.EwaldParams().DecomposedSelfEnergy( atom_self, volume );
  mprintf("DEBUG: Total self energy: %f\n", e_self);
  mprintf("DEBUG: Sum of self array: %f\n", sumArray(atom_self));
  // Do decomposed self6
  Darray const& atom_vdwself6 = NBengine_.EwaldParams().Atom_Self6Energies();
  mprintf("DEBUG: Total self6 energy: %f\n", NBengine_.EwaldParams().Self6());
  mprintf("DEBUG: Sum of self6 array: %f\n", sumArray(atom_vdwself6));

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
  mprintf("DEBUG: Recip energy      : %f\n", e_recip);
  mprintf("DEBUG: Sum of recip array: %f\n", sumArray(atom_recip));
  Darray atom_vdw6recip;
  double e_vdw6recip = LJrecip_.Recip_Decomp( atom_vdw6recip,
                                                    NBengine_.ModifyEwaldParams().SelectedCoords(),
                                                    frameIn.BoxCrd(),
                                                    NBengine_.ModifyEwaldParams().SelectedC6params(),
                                                    NBengine_.EwaldParams().LJ_EwaldCoeff()
                                                  );
  mprintf("DEBUG: VDW Recip energy      : %f\n", e_vdw6recip);
  mprintf("DEBUG: Sum of VDW recip array: %f\n", sumArray(atom_vdw6recip));
  if (NBengine_.EwaldParams().Debug() > 0) {
    mprintf("DEBUG: e_vdw6self = %16.8f\n", NBengine_.EwaldParams().Self6());
    mprintf("DEBUG: Evdwrecip = %16.8f\n", e_vdw6recip);
  }

  t_direct_.Start();
  Cpptraj::PairListTemplate<double>(pairList_, Excluded_,
                                    NBengine_.EwaldParams().Cut2(), NBengine_);
  t_direct_.Stop();
  mprintf("DEBUG: Direct Elec. energy : %f\n", NBengine_.Eelec());
  mprintf("DEBUG: Sum of elec. energy : %f\n", sumArray(NBengine_.Eatom_Elec()));
  mprintf("DEBUG: Direct VDW energy   : %f\n", NBengine_.Evdw());
  mprintf("DEBUG: Sum of VDW energy   : %f\n", sumArray(NBengine_.Eatom_EVDW()));
  mprintf("DEBUG: Direct Adjust energy: %f\n", NBengine_.Eadjust());
  mprintf("DEBUG: Sum of Adjust energy: %f\n", sumArray(NBengine_.Eatom_EAdjust()));

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
  // TODO preallocate?
  atom_elec_.resize( NBengine_.Eatom_Elec().size() );
  atom_vdw_.resize(  NBengine_.Eatom_EVDW().size() );
  for (unsigned int idx = 0; idx != atom_elec_.size(); idx++)
  {
    atom_elec_[idx] = atom_self[idx] + atom_recip[idx] + NBengine_.Eatom_Elec()[idx] + NBengine_.Eatom_EAdjust()[idx];
    atom_vdw_[idx]  = NBengine_.Eatom_EVDW()[idx] + atom_vdwself6[idx] + atom_vdw6recip[idx];
  }

  t_total_.Stop();
  return 0;
}

void EwaldCalc_Decomp_LJPME::Timing(double total) const {
  t_total_.WriteTiming(1,  "  LJPME Total:", total);
  Recip_.Timing_Total().WriteTiming(2,  "Recip:     ", t_total_.Total());
  LJrecip_.Timing_Total().WriteTiming(2,"LJRecip:   ", t_total_.Total());
  t_direct_.WriteTiming(2, "Direct:    ", t_total_.Total());
}
#endif /* LIBPME */

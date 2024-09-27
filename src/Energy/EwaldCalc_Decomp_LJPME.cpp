#include "EwaldCalc_Decomp_LJPME.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"
#include "../PairListTemplate.h"
#include "../Topology.h"

using namespace Cpptraj::Energy;

EwaldCalc_Decomp_LJPME::EwaldCalc_Decomp_LJPME() :
  Recip_(PME_Recip::COULOMB),
  LJrecip_(PME_Recip::LJ)
{}

/** Set up LJPME parameters. */
int EwaldCalc_Decomp_LJPME::Init(Box const& boxIn, EwaldOptions const& pmeOpts, int debugIn)
{
  if (NBengine_.ModifyEwaldParams().InitEwald(boxIn, pmeOpts, debugIn)) {
    mprinterr("Error: Decomposable LJPME calculation init failed.\n");
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
int EwaldCalc_Decomp_LJPME::Setup(Topology const& topIn, AtomMask const& maskIn) {
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

// DEBUG
static inline double sumArray(std::vector<double> const& arrayIn) {
  double sum = 0;
  for (std::vector<double>::const_iterator it = arrayIn.begin(); it != arrayIn.end(); ++it)
    sum += *it;
  return sum;
}

/** Calculate full nonbonded energy with LJPME */
int EwaldCalc_Decomp_LJPME::CalcDecomposedNonbondEnergy(Frame const& frameIn, AtomMask const& maskIn,
                                double& e_elec, double& e_vdw,
                                Darray& atom_elec, Darray& atom_vdw)
{
  t_total_.Start();
  double volume = frameIn.BoxCrd().CellVolume();
  Darray atom_self;
  double e_self = NBengine_.EwaldParams().DecomposedSelfEnergy( atom_self, volume );
  mprintf("DEBUG: Total self energy: %f\n", e_self);
  mprintf("DEBUG: Sum of self array: %f\n", sumArray(atom_self));
  // FIXME do decomposed self6
  Darray atom_vdwself6;

  int retVal = pairList_.CreatePairList(frameIn, frameIn.BoxCrd().UnitCell(),
                                        frameIn.BoxCrd().FracCell(), maskIn);
  if (retVal != 0) {
    mprinterr("Error: Pairlist creation failed for LJPME calc.\n");
    return 1;
  }

  // TODO make more efficient
  NBengine_.ModifyEwaldParams().FillRecipCoords( frameIn, maskIn );

  //  MapCoords(frameIn, ucell, recip, maskIn);
  Darray atom_recip;
  // FIXME helPME requires coords and charge arrays to be non-const
  double e_recip = Recip_.Recip_Decomp( atom_recip,
                                              NBengine_.ModifyEwaldParams().SelectedCoords(),
                                              frameIn.BoxCrd(),
                                              NBengine_.ModifyEwaldParams().SelectedCharges(),
                                              NBengine_.EwaldParams().NFFT(),
                                              NBengine_.EwaldParams().EwaldCoeff(),
                                              NBengine_.EwaldParams().Order()
                                            );
  mprintf("DEBUG: Recip energy      : %f\n", e_recip);
  mprintf("DEBUG: Sum of recip array: %f\n", sumArray(atom_recip));
  Darray atom_vdw6recip;
  double e_vdw6recip = LJrecip_.Recip_Decomp( atom_vdw6recip,
                                                    NBengine_.ModifyEwaldParams().SelectedCoords(),
                                                    frameIn.BoxCrd(),
                                                    NBengine_.ModifyEwaldParams().SelectedC6params(),
                                                    NBengine_.EwaldParams().NFFT(),
                                                    NBengine_.EwaldParams().LJ_EwaldCoeff(),
                                                    NBengine_.EwaldParams().Order()
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
  atom_elec.resize( NBengine_.Eatom_Elec().size() );
  atom_vdw.resize(  NBengine_.Eatom_EVDW().size() );
  for (unsigned int idx = 0; idx != atom_elec.size(); idx++)
  {
    atom_elec[idx] = atom_self[idx] + atom_recip[idx] + NBengine_.Eatom_Elec()[idx] + NBengine_.Eatom_EAdjust()[idx];
    atom_vdw[idx]  = NBengine_.Eatom_EVDW()[idx] + atom_vdwself6[idx] + atom_vdw6recip[idx];
  }

  t_total_.Stop();
  return 0;
}

void EwaldCalc_Decomp_LJPME::Timing(double total) const {
  t_total_.WriteTiming(1,  "  LJPME Total:", total);
  Recip_.Timing_Total().WriteTiming(2,  "Recip:     ", t_total_.Total());
  //Recip_.Timing_Calc().WriteTiming(3,  "Recip. Calc   :", Recip_.Timing_Total().Total());
  LJrecip_.Timing_Total().WriteTiming(2,"LJRecip:   ", t_total_.Total());
  //LJrecip_.Timing_Calc().WriteTiming(3,"LJ Recip. Calc:", LJrecip_.Timing_Total().Total());
  t_direct_.WriteTiming(2, "Direct:    ", t_total_.Total());

  pairList_.Timing(total);
}

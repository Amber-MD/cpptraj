#include "EwaldCalc_Decomp_PME.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"
#include "../Topology.h"

using namespace Cpptraj::Energy;

EwaldCalc_Decomp_PME::EwaldCalc_Decomp_PME() :
  Recip_(PME_Recip::COULOMB)
{}

/** Set up PME parameters. */
int EwaldCalc_Decomp_PME::Init(Box const& boxIn, EwaldOptions const& pmeOpts, int debugIn)
{
  if (NBengine_.ModifyEwaldParams().InitEwald(boxIn, pmeOpts, debugIn)) {
    mprinterr("Error: Decomposable PME calculation init failed.\n");
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
int EwaldCalc_Decomp_PME::Setup(Topology const& topIn, AtomMask const& maskIn) {
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
int EwaldCalc_Decomp_PME::CalcDecomposedNonbondEnergy(Frame const& frameIn, AtomMask const& maskIn,
                                     double& e_elec, double& e_vdw,
                                     Darray& atom_elec, Darray& atom_vdw)
{
  t_total_.Start();
  double volume = frameIn.BoxCrd().CellVolume();
  Darray atom_self;
  double e_self = NBengine_.EwaldParams().DecomposedSelfEnergy( atom_self, volume );
  mprintf("DEBUG: Total self energy: %f\n", e_self);

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
  mprintf("DEBUG: Recip energy: %f\n", e_recip);

  return 0;
}

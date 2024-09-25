#include "EwaldCalc_Decomp_PME.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"

using namespace Cpptraj::Energy;

EwaldCalc_Decomp_PME::EwaldCalc_Decomp_PME() :
  Recip_(PME_Recip::COULOMB)
{}

/** Set up PME parameters. */
int EwaldCalc_Decomp_PME::Init(Box const& boxIn, EwaldOptions const& pmeOpts, int debugIn)
{
  //if (NBengine_.ModifyEwaldParams().InitEwald(boxIn, pmeOpts, debugIn)) {
  //  mprinterr("Error: PME calculation init failed.\n");
  //  return 1;
  //}
  if (pairList_.InitPairList(NBengine_.EwaldParams().Cutoff(), pmeOpts.SkinNB(), debugIn))
    return 1;
  if (pairList_.SetupPairList( boxIn ))
    return 1;
  VDW_LR_.SetDebug( debugIn );
  Recip_.SetDebug( debugIn );

  return 0;
}


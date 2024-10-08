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
  if (pmeOpts.Type() == EwaldOptions::PME) {
    mprinterr("Internal Error: Options were set up for PME only.\n");
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


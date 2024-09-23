#include "Calc_LJPME.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"

using namespace Cpptraj::Energy;

Calc_LJPME::Calc_LJPME() :
  Recip_(PME_Recip::COULOMB),
  LJrecip_(PME_Recip::LJ)
{}

/** Set up LJPME parameters. */
int Calc_LJPME::Init(Box const& boxIn, EwaldOptions const& pmeOpts, int debugIn)
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


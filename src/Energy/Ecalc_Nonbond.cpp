#include "Ecalc_Nonbond.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
Ecalc_Nonbond::Ecalc_Nonbond() :
  type_(UNSPECIFIED),
  needs_pairlist_(false)
{}

/** Init */
int Ecalc_Nonbond::InitNonbondCalc(CalcType typeIn, Box const& boxIn,
                                   EwaldOptions const& pmeOpts, int debugIn)
{
  if (typeIn == UNSPECIFIED) {
    mprinterr("Internal Error: Ecalc_Nonbond::InitNonbondCalc(): No nonbonded calc type specified.\n");
    return 1;
  }
  needs_pairlist_ = false;

  if (needs_pairlist_) {
    if (pairList_.InitPairList(pmeOpts.Cutoff(), pmeOpts.SkinNB(), debugIn))
      return 1;
    if (pairList_.SetupPairList( boxIn ))
      return 1;
  }

  return 0;
}

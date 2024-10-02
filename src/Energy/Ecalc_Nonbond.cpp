#include "Ecalc_Nonbond.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"
#include "../Topology.h"

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

/** Setup */
int Ecalc_Nonbond::SetupNonbondCalc(Topology const& topIn, AtomMask const& maskIn) {
  // Setup exclusion list
  // Use distance of 4 (up to dihedrals)
  if (Excluded_.SetupExcluded(topIn.Atoms(), maskIn, 4,
                              ExclusionArray::EXCLUDE_SELF,
                              ExclusionArray::FULL))
  {
    mprinterr("Error: Could not set up exclusion list for nonbonded calculation.\n");
    return 1;
  }
  return 0;
}

/** Calc */
int Ecalc_Nonbond::NonbondEnergy(Frame const& frameIn, AtomMask const& maskIn,
                                 double& e_elec, double& e_vdw)
{
  if (needs_pairlist_) {
    if (pairList_.CreatePairList(frameIn, frameIn.BoxCrd().UnitCell(),
                                 frameIn.BoxCrd().FracCell(), maskIn) != 0)
    {
      mprinterr("Error: Pairlist creation failed for nonbond calc.\n");
      return 1;
    }
  }

  return 0;
}

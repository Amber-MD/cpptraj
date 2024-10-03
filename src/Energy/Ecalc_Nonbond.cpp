#include "Ecalc_Nonbond.h"
#include "EwaldCalc_PME.h"
#include "EwaldCalc_LJPME.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"
#include "../Topology.h"
#include <cmath> // sqrt for Ene_Nonbond
#include "Ene_Nonbond.h"

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
Ecalc_Nonbond::Ecalc_Nonbond() :
  calc_(0),
  currentTop_(0),
  type_(UNSPECIFIED),
  needs_pairlist_(false)
{}

/** DESTRUCTOR */
Ecalc_Nonbond::~Ecalc_Nonbond() {
  if (calc_ != 0)
    delete calc_;
}

/** Init */
int Ecalc_Nonbond::InitNonbondCalc(CalcType typeIn, Box const& boxIn,
                                   EwaldOptions const& pmeOpts, int debugIn)
{
  currentTop_ = 0;
  type_ = typeIn;

  needs_pairlist_ = false;
  calc_ = 0;
  switch (type_) {
    case SIMPLE :
      break;
    case PME    :
      needs_pairlist_ = true;
      calc_ = new EwaldCalc_PME();
      break;
    case LJPME  :
      needs_pairlist_ = true;
      calc_ = new EwaldCalc_LJPME();
      break;
    case UNSPECIFIED :
      mprinterr("Internal Error: Ecalc_Nonbond::InitNonbondCalc(): No nonbonded calc type specified.\n");
      return 1;
  }

  if (needs_pairlist_) {
    if (pairList_.InitPairList(pmeOpts.Cutoff(), pmeOpts.SkinNB(), debugIn))
      return 1;
    if (pairList_.SetupPairList( boxIn ))
      return 1;
  }
  if (calc_ != 0) {
    if (calc_->Init( boxIn, pmeOpts, debugIn ))
      return 1;
  }

  return 0;
}

/** Setup */
int Ecalc_Nonbond::SetupNonbondCalc(Topology const& topIn, AtomMask const& maskIn) {
  currentTop_ = &topIn;
  // Setup exclusion list
  // Use distance of 4 (up to dihedrals)
  if (Excluded_.SetupExcluded(topIn.Atoms(), maskIn, 4,
                              ExclusionArray::EXCLUDE_SELF,
                              ExclusionArray::FULL))
  {
    mprinterr("Error: Could not set up exclusion list for nonbonded calculation.\n");
    return 1;
  }
  if (calc_ != 0) {
    if (calc_->Setup(topIn, maskIn))
      return 1;
  }
  return 0;
}

/** Calc */
int Ecalc_Nonbond::NonbondEnergy(Frame const& frameIn, AtomMask const& maskIn,
                                 double& e_elec, double& e_vdw)
{
  t_total_.Start();

  int err = 0;
  if (type_ == SIMPLE) {
    static const double QFAC = Constants::ELECTOAMBER * Constants::ELECTOAMBER;
    Ene_Nonbond<double>(frameIn, *currentTop_, maskIn, Excluded_, QFAC,
                        e_elec, e_vdw);
  } else {

    if (needs_pairlist_) {
      if (pairList_.CreatePairList(frameIn, frameIn.BoxCrd().UnitCell(),
                                   frameIn.BoxCrd().FracCell(), maskIn) != 0)
      {
        mprinterr("Error: Pairlist creation failed for nonbond calc.\n");
        return 1;
      }
    }

    err = calc_->CalcNonbondEnergy(frameIn, maskIn, pairList_, Excluded_,
                                   e_elec, e_vdw);
  }
  t_total_.Stop();
  return err;
}

/** Print timing */
void Ecalc_Nonbond::PrintTiming(double total) const {
  if (calc_ != 0) calc_->Timing( t_total_.Total() );
  if (needs_pairlist_)
    pairList_.Timing( t_total_.Total() );
  t_total_.WriteTiming(0, "Nonbond total:");
}

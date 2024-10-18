#include "Ecalc_Nonbond.h"
#include "EwaldCalc_Regular.h"
#ifdef LIBPME
# include "EwaldCalc_LJPME.h"
# include "EwaldCalc_PME.h"
# include "EwaldCalc_Decomp_LJPME.h"
# include "EwaldCalc_Decomp_PME.h"
#else
// TODO This is here for when LIBPME is not defined because there is
//      currently no decomp version of regular Ewald. If there ever is
//      this inlclude can be removed.
# include "EwaldCalc_Decomp.h"
#endif
#include "../CharMask.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"
#include "../Topology.h"
#include <cmath> // sqrt for Ene_Nonbond, Ene_Decomp_Nonbond
#include "Ene_Nonbond.h"
#include "Ene_Decomp_Nonbond.h"

using namespace Cpptraj::Energy;

/** CONSTRUCTOR */
Ecalc_Nonbond::Ecalc_Nonbond() :
  calc_(0),
  currentTop_(0),
  type_(UNSPECIFIED),
  decompose_energy_(false)
{}

/** DESTRUCTOR */
Ecalc_Nonbond::~Ecalc_Nonbond() {
  if (calc_ != 0)
    delete calc_;
}

/** Init */
int Ecalc_Nonbond::InitNonbondCalc(CalcType typeIn, bool decompose_energyIn,
                                   Box const& boxIn, EwaldOptions const& pmeOpts, int debugIn)
{
  currentTop_ = 0;
  type_ = typeIn;
  decompose_energy_ = decompose_energyIn;

  calc_ = 0;
  switch (type_) {
    case SIMPLE :
      break;
#   ifdef LIBPME
    case PME    :
      if (decompose_energy_)
        calc_ = new EwaldCalc_Decomp_PME();
      else
        calc_ = new EwaldCalc_PME();
      break;
    case LJPME  :
      if (decompose_energy_)
        calc_ = new EwaldCalc_Decomp_LJPME();
      else
        calc_ = new EwaldCalc_LJPME();
      break;
#   else
    case PME : mprinterr("Error: PME requires compiling with LIBPME (FFTW3 and C++11 support).\n"); return 1;
    case LJPME : mprinterr("Error: LJPME requires compiling with LIBPME (FFTW3 and C++11 support).\n"); return 1;
#   endif
    case REGULAR_EWALD :
      if (decompose_energy_) {
        mprinterr("Internal Error: Ecalc_Nonbond::InitNonbondCalc(): Cannot decompose regular Ewald calc.\n");
        return 1;
      } else
        calc_ = new EwaldCalc_Regular();
      break;
    case UNSPECIFIED :
      mprinterr("Internal Error: Ecalc_Nonbond::InitNonbondCalc(): No nonbonded calc type specified.\n");
      return 1;
  }
  if (type_ != SIMPLE && calc_ == 0) {
    mprinterr("Internal Error: Ecalc_Nonbond::InitNonbondCalc(): Ewald calc alloc failed.\n");
    return 1;
  }

  if (type_ != SIMPLE && !boxIn.HasBox()) {
    mprinterr("Error: Ewald requires unit cell information.\n");
    return 1;
  }

  if (calc_ != 0) {
    // Ewald calcs need pairlist
    if (pairList_.InitPairList(pmeOpts.Cutoff(), pmeOpts.SkinNB(), debugIn))
      return 1;
    if (pairList_.SetupPairList( boxIn ))
      return 1;
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
    if (calc_->Setup(topIn, maskIn)) {
      mprinterr("Error: Nonbonded calculation setup failed.\n");
      return 1;
    }
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
    Ene_Nonbond<double>(frameIn, *currentTop_, maskIn, Excluded_, Constants::COULOMBFACTOR,
                        e_elec, e_vdw);
  } else if (calc_ != 0) {
    if (pairList_.CreatePairList(frameIn, frameIn.BoxCrd().UnitCell(),
                                 frameIn.BoxCrd().FracCell(), maskIn) != 0)
    {
      mprinterr("Error: Pairlist creation failed for nonbond calc.\n");
      return 1;
    }

    err = calc_->CalcNonbondEnergy(frameIn, maskIn, pairList_, Excluded_,
                                   e_elec, e_vdw);
  } else
    return 1; // Sanity check

  t_total_.Stop();
  return err;
}

/** Energy decomp calc */
int Ecalc_Nonbond::DecomposedNonbondEnergy(Frame const& frameIn, CharMask const& cmaskIn,
                                           double& e_elec, double& e_vdw,
                                           Darray& atom_elec, Darray& atom_vdw)
{
  // sanity check
  if (!decompose_energy_) {
    mprinterr("Internal Error: Ecalc_Nonbond::DecomposedNonbondEnergy(): Not set up for energy decomp.\n");
    return 1;
  }
  t_total_.Start();

  int err = 0;
  if (type_ == SIMPLE) {
    Ene_Decomp_Nonbond<double>(frameIn, *currentTop_, cmaskIn, Excluded_, Constants::COULOMBFACTOR,
                               e_elec, e_vdw, atom_elec, atom_vdw);
  } else if (calc_ != 0) {
    // FIXME this is an unneeded atom mask.
    AtomMask tmpMask(0, frameIn.Natom());

    if (pairList_.CreatePairList(frameIn, frameIn.BoxCrd().UnitCell(),
                                 frameIn.BoxCrd().FracCell(), tmpMask) != 0)
    {
      mprinterr("Error: Pairlist creation failed for nonbond calc.\n");
      return 1;
    }

    err = calc_->CalcNonbondEnergy(frameIn, tmpMask, pairList_, Excluded_,
                                   e_elec, e_vdw);
    EwaldCalc_Decomp const* EW_ENE = static_cast<EwaldCalc_Decomp const*>( calc_ );
    for (int at = 0; at < frameIn.Natom(); at++) {
      if (cmaskIn.AtomInCharMask(at)) {
        atom_elec[at] += EW_ENE->Atom_Elec()[at];
        atom_vdw[at] += EW_ENE->Atom_VDW()[at];
      }
    }
  } else
    return 1; // sanity check

  t_total_.Stop();
  return err;
}


/** Print timing */
void Ecalc_Nonbond::PrintTiming(double total) const {
  t_total_.WriteTiming(1,    "NONBOND     :", total);
  if (calc_ != 0) {
    calc_->Timing( t_total_.Total() );
    pairList_.Timing( t_total_.Total() );
  }
}

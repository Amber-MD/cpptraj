#include "PotentialFunction.h"
#include "CpptrajStdio.h"
#include "Topology.h"
#include "MdOpts.h"
// ----- All potential terms -----------
#include "PotentialTerm_Bond.h"
#include "PotentialTerm_OpenMM.h"
#include "PotentialTerm_LJ_Coulomb.h"
#include "PotentialTerm_Angle.h"
#include "PotentialTerm_Dihedral.h"

/** Add a term to the potential function. */
int PotentialFunction::AddTerm(PotentialTerm::Type typeIn, MdOpts const& opts) {
  PotentialTerm* term = 0;
  switch (typeIn) {
    case PotentialTerm::BOND : term = (PotentialTerm*)new PotentialTerm_Bond(); break;
    case PotentialTerm::OPENMM : term = (PotentialTerm*)new PotentialTerm_OpenMM(); break;
    case PotentialTerm::SIMPLE_LJ_Q : term = (PotentialTerm*)new PotentialTerm_LJ_Coulomb(); break;
    case PotentialTerm::ANGLE : term = (PotentialTerm*)new PotentialTerm_Angle(); break;
    case PotentialTerm::DIHEDRAL : term = (PotentialTerm*)new PotentialTerm_Dihedral(); break;
    default :
      mprinterr("Internal Error: No allocator type for potential term type #%i.\n", (int)typeIn);
      return 1;
  }
  if (term == 0) {
    mprinterr("Internal Error: Could not allocate potential term %s.\n", PotentialTerm::TypeStr(typeIn));
    return 1;
  }
  if (term->InitTerm(opts)) {
    mprinterr("Error: Term init failed for %s.\n", term->TypeStr());
    return 1;
  }
  terms_.push_back( term );
  return 0;
}

/** Initialize (or re-initialize) each term with given options. */
int PotentialFunction::InitPotential(MdOpts const& optsIn) {
  for (Parray::const_iterator it = terms_.begin(); it != terms_.end(); ++it)
  {
    if ( (*it)->InitTerm( optsIn ) ) {
      mprinterr("Error: Could not initialize term %s\n", (*it)->TypeStr());
      return 1;
    }
  }
  return 0;
}

/** Set up each term of the potential function. */
int PotentialFunction::SetupPotential(Topology const& topIn, Box const& boxIn,
                                      std::string const& maskExprIn)
{
  // First set up the mask
  mask_.ResetMask();
  if (mask_.SetMaskString( maskExprIn )) {
    mprinterr("Error: Could not set up mask expression.\n");
    return 1;
  }
  if (topIn.SetupCharMask( mask_ )) {
    mprinterr("Error: Could not set up mask.\n");
    return 1;
  }
  mask_.MaskInfo();
  return setupPotential(topIn, boxIn);
}

/** Add term with default options. */
int PotentialFunction::AddTerm(PotentialTerm::Type typeIn) {
  MdOpts opts;
  mprintf("\tUsing default options for term %s.\n", PotentialTerm::TypeStr(typeIn));
  return AddTerm(typeIn, opts);
}

/** Set up potential. */
int PotentialFunction::SetupPotential(Topology const& topIn, Box const& boxIn, CharMask const& maskIn) {
  if (maskIn.Nselected() < 1) {
    mprinterr("Internal Error: SetupPotential called with empty mask.\n");
    return 1;
  }
  mask_ = maskIn;
  return setupPotential(topIn, boxIn);
}

int PotentialFunction::setupPotential(Topology const& topIn, Box const& boxIn) {
  // Determine degrees of freedom
  deg_of_freedom_ = 3 * mask_.Nselected();
  mprintf("\t%i degrees of freedom.\n", deg_of_freedom_);

  earray_.clear();
  for (Parray::const_iterator it = terms_.begin(); it != terms_.end(); ++it)
  {
    //mprintf("DEBUG: Set up term %s\n", (*it)->TypeStr());
    if ( (*it)->SetupTerm( topIn, boxIn, mask_, earray_ ) ) {
      mprinterr("Error: Could not set up potential function term %s.\n", (*it)->TypeStr());
      return 1;
    }
    deg_of_freedom_ -= (*it)->RemovedDegreesOfFreedom();
  }
  mprintf("\tCorrected degrees of freedom: %i\n", deg_of_freedom_);
  current_ = (Topology*)&(topIn);
  return 0;
}

/** Calculate force for each term. */
int PotentialFunction::CalculateForce(Frame& frameIn) {
  earray_.zero();
  frameIn.ZeroForces();
  for (Parray::const_iterator it = terms_.begin(); it != terms_.end(); ++it)
  {
    (*it)->CalcForce( frameIn, mask_ );
  }
  return 0;
}

/** Print information about potential function. */
void PotentialFunction::FnInfo() const {
  mprintf("\t%zu terms:", terms_.size());
  for (Parray::const_iterator it = terms_.begin(); it != terms_.end(); ++it)
    mprintf(" %s", (*it)->TypeStr());
  mprintf("\n");
  mprintf("\tEnergy terms:");
  earray_.PrintActiveTerms();
  mprintf("\n");
}

#include "PotentialFunction.h"
#include "CpptrajStdio.h"
#include "Topology.h"
#include "Box.h"
// ----- All potential terms -----------
#include "PotentialTerm_Bond.h"

/** Add a term to the potential function. */
int PotentialFunction::AddTerm(PotentialTerm::Type typeIn) {
  PotentialTerm* term = 0;
  switch (typeIn) {
    case PotentialTerm::BOND : term = (PotentialTerm*)new PotentialTerm_Bond(); break;
    default :
      mprinterr("Internal Error: No allocator type for potential term.\n");
      return 1;
  }
  if (term == 0) {
    mprinterr("Internal Error: Could not allocate potential term.\n");
    return 1;
  }
  terms_.push_back( term );
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

  // Determine degrees of freedom
  // TODO depending on what terms are present and how they are set up the DoF calc may change
  deg_of_freedom_ = 3 * mask_.Nselected();
  mprintf("\t%i degrees of freedom.\n", deg_of_freedom_);

  earray_.clear();
  for (Parray::const_iterator it = terms_.begin(); it != terms_.end(); ++it)
  {
    if ( (*it)->SetupTerm( topIn, mask_, earray_ ) ) {
      mprinterr("Error: Could not set up energy term.\n");
      return 1;
    }
  }
  current_ = (Topology*)&(topIn);
  return 0;
}

/** Calculate force for each term. */
int PotentialFunction::CalculateForce(Frame& frameIn) {
  earray_.zero();
  for (Parray::const_iterator it = terms_.begin(); it != terms_.end(); ++it)
  {
    (*it)->CalcForce( frameIn, mask_ );
  }
  return 0;
}

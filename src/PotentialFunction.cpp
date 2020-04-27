#include "PotentialFunction.h"
#include "CpptrajStdio.h"
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
int PotentialFunction::SetupPotential(Topology const& topIn, CharMask const& maskIn) {
  earray_.clear();
  for (Parray::const_iterator it = terms_.begin(); it != terms_.end(); ++it)
  {
    if ( (*it)->SetupTerm( topIn, maskIn, earray_ ) ) {
      mprinterr("Error: Could not set up energy term.\n");
      return 1;
    }
  }
  return 0;
}

/** Calculate force for each term. */
int PotentialFunction::CalculateForce(Frame& frameIn, CharMask const& maskIn) {
  earray_.zero();
  for (Parray::const_iterator it = terms_.begin(); it != terms_.end(); ++it)
  {
    (*it)->CalcForce( frameIn, maskIn );
  }
  return 0;
}

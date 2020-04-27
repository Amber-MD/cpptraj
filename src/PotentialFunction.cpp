#include "PotentialFunction.h"
#include "CpptrajStdio.h"
// ----- All potential terms -----------
#include "PotentialTerm_Bond.h"

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

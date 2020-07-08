#include "PotentialTerm_InitOpts.h"

PotentialTerm::InitOpts::InitOpts() :
  constraintType_(NO_CONSTRAINTS),
  scaleEE_(1.0/1.2), // Amber default
  scaleNB_(1.0/2.0), // Amber default
  cutEE_(8.0),       // in Ang., Amber default
  cutNB_(8.0)        // in Ang., Amber default
{}

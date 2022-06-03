#include "FunctionalGroup.h"

using namespace Cpptraj::Structure;

FunctionalGroup::FunctionalGroup() {}

/** Corresponds to FunctionalGroupType */
const char* FunctionalGroup::FunctionalGroupStr_[] = {
  "SO3", "CH3", "Acetyl", "OH", "OCH3", "Unrecognized"
};


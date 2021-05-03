#include "Pucker_PuckerMask.h"

using namespace Cpptraj;

Pucker::PuckerMask::PuckerMask() {}

Pucker::PuckerMask::PuckerMask(std::vector<int> const& atomsIn) :
  atoms_(atomsIn)
{}

#include "Pucker_PuckerMask.h"
#include "Topology.h"

using namespace Cpptraj;

Pucker::PuckerMask::PuckerMask() {}

Pucker::PuckerMask::PuckerMask(std::vector<int> const& atomsIn) :
  atoms_(atomsIn)
{}

std::string Pucker::PuckerMask::PuckerMaskString(Topology const& topIn) const {
  std::string out;
  for (std::vector<int>::const_iterator it = atoms_.begin(); it != atoms_.end(); ++it)
  {
    if (it != atoms_.begin())
      out.append(" ");
    out.append( topIn.TruncResNameAtomName( *it ) );
  }
  return out;
}

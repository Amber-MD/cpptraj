#include "Pucker_PuckerMask.h"
#include "Topology.h"

using namespace Cpptraj;

Pucker::PuckerMask::PuckerMask() {}

/** CONSTRUCTOR - res #, name, index */
Pucker::PuckerMask::PuckerMask(int resnumIn, std::string const& nameIn, std::vector<int> const& atomsIn) :
  atoms_(atomsIn),
  aspect_(nameIn),
  resnum_(resnumIn)
{}

/** \return string based on atoms in mask. */
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

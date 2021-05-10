#include "Pucker_PuckerToken.h"
#include "Topology.h"
#include "Pucker_PuckerMask.h"

using namespace Cpptraj;

Pucker::PuckerToken::PuckerToken() {}

Pucker::PuckerToken::PuckerToken(std::string const& nameIn, NameArray const& anamesIn) :
  name_(nameIn), atomNames_(anamesIn)
{}

/** Recursive function for finding pucker atoms. */
void Pucker::PuckerToken::FindAtoms(Topology const& topIn, int at,
                                    unsigned int idx, unsigned int maxidx, std::vector<int>& indices)
const
{
  if (idx >= maxidx) return;

  for (Atom::bond_iterator it = topIn[at].bondbegin(); it != topIn[at].bondend(); ++it)
  {
    if (topIn[*it].Name() == atomNames_[idx]) {
      indices[idx] = *it;

      FindAtoms( topIn, *it, idx+1, maxidx, indices );
    }
  }
}

/** \return PuckerMask containing indices of this pucker found in specified residue.
  */
Pucker::PuckerMask Pucker::PuckerToken::FindPuckerAtoms(Topology const& topIn, int resnum)
const
{
  Residue const& currentRes = topIn.Res( resnum );

  std::vector<int> indices(atomNames_.size(), -1);

  // Find the first atom
  for (int at = currentRes.FirstAtom(); at != currentRes.LastAtom(); ++at)
  {
    if (topIn[at].Name() == atomNames_.front()) {
      indices[0] = at;

      FindAtoms( topIn, at, 1, atomNames_.size(), indices );
      break;
    }
  }
  return PuckerMask( indices );
}

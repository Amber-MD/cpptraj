#include "Pucker_PuckerToken.h"
#include "Topology.h"
#include "Pucker_PuckerMask.h"
#include "CpptrajStdio.h"

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
  if (idx >= maxidx) {
    // Make sure this is a cycle (i.e. this should be bonded to the first atom detected)
    //mprintf("DEBUG: Checking for cycle at=%i idx=%u maxidx=%u indices=", at, idx, maxidx);
    //for (std::vector<int>::const_iterator it = indices.begin(); it != indices.end(); ++it)
    //  mprintf(" %i", *it);
    //mprintf("\n");
    if (idx > maxidx) {
      // Sanity check
      mprinterr("Internal Error: PuckerToken::FindAtoms: Index out of range.\n");
      return;
    }
    for (Atom::bond_iterator it = topIn[at].bondbegin(); it != topIn[at].bondend(); ++it)
    {
      if (*it == indices[0]) {
        indices[idx] = *it;
        return;
      }
    }
  } else {
    // Check that this atom is bonded to the next one expected in the sequence.
    for (Atom::bond_iterator it = topIn[at].bondbegin(); it != topIn[at].bondend(); ++it)
    {
      if (topIn[*it].Name() == atomNames_[idx]) {
        indices[idx] = *it;

        FindAtoms( topIn, *it, idx+1, maxidx, indices );
      }
    }
  }
}

/** \return PuckerMask containing indices of this pucker found in specified residue.
  */
Pucker::PuckerMask Pucker::PuckerToken::FindPuckerAtoms(Topology const& topIn, int resnum)
const
{
  Residue const& currentRes = topIn.Res( resnum );

  std::vector<int> indices(atomNames_.size()+1, -1);

  // Find the first atom
  for (int at = currentRes.FirstAtom(); at != currentRes.LastAtom(); ++at)
  {
    if (topIn[at].Name() == atomNames_.front()) {
      indices[0] = at;

      FindAtoms( topIn, at, 1, atomNames_.size(), indices );
      break;
    }
  }
  // DEBUG
  //mprintf("DEBUG: Results for pucker '%s', residue %i:", name_.c_str(), resnum+1);
  //for (std::vector<int>::const_iterator it = indices.begin(); it != indices.end(); ++it)
  //  mprintf(" %i", *it);
  //mprintf("\n");
  // Check if the entire pucker was found.
  std::vector<int> actualIndices;
  actualIndices.reserve(atomNames_.size());
  if (indices.back() == -1) {
    // This means a cycle was not found.
    return PuckerMask();
  }
  for (unsigned int idx = 0; idx != atomNames_.size(); idx++)
  {
    if (indices[idx] == -1) {
      return PuckerMask();
    }
    actualIndices.push_back( indices[idx] );
  }
  return PuckerMask( resnum, Name(), actualIndices );
}

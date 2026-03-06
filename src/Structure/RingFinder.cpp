#include "RingFinder.h"
#include "../AtomMask.h"
#include "../CharMask.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"
#include <map>

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
RingFinder::RingFinder() :
  debug_(0)
{}

/** Set debug level */
void RingFinder::SetDebug(int debugIn) {
  debug_ = debugIn;
}

/** Init ring finder. */
int RingFinder::InitRingFinder(ArgList& argIn) {
  return 0;
}

static void visitAtom( int at, int previousAt, int res, Topology const& topIn, std::vector<bool>& Visited, std::vector<int>& startAtoms )
{
  if (Visited[at]) {
    //mprintf("Already visited %s %zu\n", topIn.AtomMaskName(at).c_str(), startAtoms.size());
    startAtoms.push_back( at );
  } else {
    Visited[at] = true;
    Atom const& currentAt = topIn[at];
    for (Atom::bond_iterator bat = currentAt.bondbegin(); bat != currentAt.bondend(); ++bat)
    {
      if (*bat != previousAt && topIn[*bat].Nbonds() > 1 && topIn[*bat].ResNum() == res) {
        visitAtom( *bat, at, res, topIn, Visited, startAtoms );
      }
    }
  }
}

static std::vector<int> findCycle(int at, int tgtAt, int previousAt, int res, Topology const& topIn,
                                  std::vector<bool> Visited, std::vector<int> ringAtoms)
{
  if (at == tgtAt) {
    //mprintf("CYCLE FOUND (%zu).", ringAtoms.size());
    //for (std::vector<int>::const_iterator it = ringAtoms.begin(); it != ringAtoms.end(); ++it)
    //  mprintf(" %s", topIn.AtomMaskName(*it).c_str());
    //mprintf("\n");
    return ringAtoms;
  }
  Visited[at] = true;
  ringAtoms.push_back( at );
  Atom const& currentAt = topIn[at];
  typedef std::vector<int> Iarray;
  std::vector<Iarray> ringAtArray( currentAt.Nbonds() );
  std::vector<Iarray>::iterator ringIt = ringAtArray.begin();
  std::vector<Iarray>::iterator shortestRing = ringAtArray.end();
  for (Atom::bond_iterator bat = currentAt.bondbegin(); bat != currentAt.bondend(); ++bat, ++ringIt)
  {
    if (*bat != previousAt && topIn[*bat].Nbonds() > 1 &&
        topIn[*bat].ResNum() == res && !Visited[*bat])
    {
      *ringIt = findCycle(*bat, tgtAt, at, res, topIn, Visited, ringAtoms);
      if (!ringIt->empty()) {
        if (shortestRing == ringAtArray.end())
          shortestRing = ringIt;
        else if (ringIt->size() < shortestRing->size())
          shortestRing = ringIt;
      }
    }
  }
  // If multiple paths were found, return the shortest.
  if (shortestRing != ringAtArray.end())
    return *shortestRing;
  return std::vector<int>();
}

void RingFinder::useCachedRings(Topology const& topIn, int res, CharMask const& cmask, RingNamesType const& ringNamesArray)
{
  if (debug_ > 0)
    mprintf("DEBUG: Using cached values for residue %s, %zu rings.\n", topIn.TruncResNameNum(res).c_str(), ringNamesArray.size());
  for (RingNamesType::const_iterator rnames = ringNamesArray.begin(); rnames != ringNamesArray.end(); ++rnames)
  {
    // Find all ring atoms. They must be present and selected.
    bool allAtomsFound = true;
    std::vector<int> ringAtoms;
    ringAtoms.reserve(rnames->size());
    for (AtomNameArray::const_iterator aname = rnames->begin(); aname != rnames->end(); ++aname)
    {
      NameType const& atmName = *aname;
      int at = topIn.FindAtomInResidue(res, atmName);
      if (at < 0) {
        mprintf("Warning: Ring atom %s not found in residue %s\n", *atmName, topIn.TruncResNameNum(res).c_str());
        allAtomsFound = false;
        break;
      }
      if (!cmask.AtomInCharMask(at)) {
        allAtomsFound = false;
        break;
      }
      ringAtoms.push_back(at);
    }
    if (allAtomsFound)
      rings_.push_back( AtomMask(ringAtoms, topIn.Natom()) );
  } // END loop over rings in residue
}

/** Set up ring finder for topology. */
int RingFinder::SetupRingFinder(Topology const& topIn, AtomMask const& maskIn) {
  rings_.clear();

  typedef std::pair<NameType, RingNamesType> RingCachePair;
  typedef std::map<NameType, RingNamesType> RingCacheType;
  // Hold names of ring atoms for each unique type of residue.
  RingCacheType RingCache;

  CharMask cmask( maskIn.ConvertToCharMask(), maskIn.Nselected() );

  std::vector<bool> Visited(topIn.Natom(), false); // TODO only have visited for residue atoms
  // Do not look for rings that span residues.
  for (int res = 0; res < topIn.Nres(); res++) {
    Residue const& currentRes = topIn.Res(res);
    // Check the cache
    RingCacheType::const_iterator it = RingCache.find( currentRes.Name() );
    if ( it != RingCache.end() ) {
      useCachedRings( topIn, res, cmask, it->second );
      continue;
    }
    // First time this residue has been seen. Do the full ring search.
    RingNamesType ringNamesArray;
    // Ignore solvent
    if (topIn.Mol( topIn[currentRes.FirstAtom()].MolNum() ).IsSolvent()) continue;
    std::vector<int> startAtoms;
    int idx = 0;
    // Initial scan for start atoms
    for (int at = currentRes.FirstAtom(); at != currentRes.LastAtom(); at++, idx++)
    {
      if (cmask.AtomInCharMask( at )) {
        Atom const& currentAt = topIn[at];
        if (currentAt.Nbonds() > 1) {
          //mprintf("Starting at %s\n", topIn.AtomMaskName(at).c_str());
          Visited.assign( topIn.Natom(), false );
          Visited[at] = true;
          for (Atom::bond_iterator bat = currentAt.bondbegin(); bat != currentAt.bondend(); ++bat)
          {
            if (topIn[*bat].Nbonds() > 1 && topIn[*bat].ResNum() == res) {
              visitAtom( *bat, at, res, topIn, Visited, startAtoms );
              if (!startAtoms.empty()) break;
            }
          }
        }
        if (!startAtoms.empty()) break;
      }
    }
    if (debug_ > 0 && !startAtoms.empty()) {
      mprintf("Potential start atoms:");
      for (std::vector<int>::const_iterator it = startAtoms.begin(); it != startAtoms.end(); ++it)
        mprintf(" %s", topIn.AtomMaskName(*it).c_str());
      mprintf("\n");
    }
    // Loop over potential start atoms
    Visited.assign( topIn.Natom(), false );
    for (std::vector<int>::const_iterator it = startAtoms.begin(); it != startAtoms.end(); ++it)
    {
      if (!cmask.AtomInCharMask( *it )) continue;
      Atom const& startAt = topIn[*it];
      // Need to have at least 2 bonds
      if (startAt.Nbonds() > 1) {
        // Loop over all other combinations of bonded atoms. Want to see if
        // we can make a cycle from at1 to at2:
        // at1 - startAt - at2 TODO check that bonded atoms are in the same residue
        for (int bidx1 = 0; bidx1 != startAt.Nbonds(); bidx1++) {
          int at1 = startAt.Bond(bidx1);
          if (topIn[at1].ResNum() == res) {
            for (int bidx2 = bidx1 + 1; bidx2 != startAt.Nbonds(); bidx2++) {
              int at2 = startAt.Bond(bidx2);
              if (topIn[at2].ResNum() == res) {
                //mprintf("Check <- %s - %s - %s ->\n",
                //        topIn.AtomMaskName(at1).c_str(),
                //        topIn.AtomMaskName(*it).c_str(),
                //        topIn.AtomMaskName(at2).c_str());
                Visited[*it] = true;
                std::vector<int> tmpAtoms;
                std::vector<int> ringAtoms = findCycle(at1, at2, *it, res, topIn, Visited, tmpAtoms);
                if (!ringAtoms.empty()) {
                  // Fill in the rest of the atoms
                  ringAtoms.push_back( at2 );
                  ringAtoms.push_back( *it );
                  // Mark the ring atoms as visited
                  Visited[at2] = true;
                  for (std::vector<int>::const_iterator rt = ringAtoms.begin(); rt != ringAtoms.end(); ++rt)
                    Visited[*rt] = true;
                  // Only add the ring if all atoms were found
                  bool allfound = true;
                  if (debug_ > 0) mprintf("RING FOUND (%zu): {", ringAtoms.size());
                  for (std::vector<int>::const_iterator rt = ringAtoms.begin(); rt != ringAtoms.end(); ++rt) {
                    if (debug_ > 0) mprintf(" %s", topIn.AtomMaskName(*rt).c_str());
                    if (!cmask.AtomInCharMask( *rt )) allfound = false;
                  }
                  if (debug_ > 0) mprintf(" }\n");
                  if (allfound) {
                    rings_.push_back( AtomMask(ringAtoms, topIn.Natom()) );
                    AtomNameArray rnames;
                    rnames.reserve( ringAtoms.size() );
                    for (std::vector<int>::const_iterator zt = ringAtoms.begin(); zt != ringAtoms.end(); ++zt)
                      rnames.push_back( topIn[*zt].Name() );
                    ringNamesArray.push_back( rnames );
                  }
                }
              } // END bond atom 2 in residue
            } // END inner loop over bonds
          } // END bond atom 1 in residue
        } // END outer loop over bonds
      } // END # bonds > 1
    } // END loop over potential start atoms
    RingCache.insert( RingCachePair( currentRes.Name(), ringNamesArray ) );
  } // END loop over all residues

  return 0;
}

/** Print found rings to stdout. */
void RingFinder::PrintRings(Topology const& topIn) const {
  mprintf("\t%zu rings found.\n", rings_.size());
  for (Marray::const_iterator mask = rings_.begin(); mask != rings_.end(); ++mask) {
    mprintf("\t  %i atoms:", mask->Nselected());
    for (AtomMask::const_iterator it = mask->begin(); it != mask->end(); ++it)
      mprintf(" %s", topIn.AtomMaskName(*it).c_str());
    mprintf("\n");
  }
}

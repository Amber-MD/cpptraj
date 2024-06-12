#include "HisProt.h"
#include "StructureRoutines.h"
#include "../CpptrajStdio.h"
#include "../NameType.h"
#include "../Topology.h"

using namespace Cpptraj::Structure;

/** Try to determine histidine protonation from existing hydrogens.
  * Change residue names as appropriate.
  * \param topIn Input topology.
  * \param ND1 ND1 atom name.
  * \param NE2 NE2 atom name.
  * \param HisName PDB histidine name.
  * \param HieName Name for epsilon-protonated His.
  * \param HidName Name for delta-protonated His.
  * \param HipName Name for doubly-protonated His.
  */
int Cpptraj::Structure::DetermineHisProt( Topology& topIn,
                                           NameType const& ND1,
                                           NameType const& NE2,
                                           NameType const& HisName,
                                           NameType const& HieName,
                                           NameType const& HidName,
                                           NameType const& HipName )
{
  typedef std::vector<int> Iarray;
  mprintf("\tAttempting to determine histidine form from any existing H atoms.\n");
  std::string hisMaskStr = ":" + HisName.Truncated();
  AtomMask mask;
  if (mask.SetMaskString( hisMaskStr )) {
    mprinterr("Error: Invalid His mask string: %s\n", hisMaskStr.c_str());
    return 1;
  }
  if ( topIn.SetupIntegerMask( mask )) return 1;
  mask.MaskInfo();
  Iarray resIdxs = topIn.ResnumsSelectedBy( mask );
  // Loop over selected histidine residues
  unsigned int nchanged = 0;
  for (Iarray::const_iterator rnum = resIdxs.begin();
                              rnum != resIdxs.end(); ++rnum)
  {
    if (StructureDebugLevel() > 1)
      mprintf("DEBUG: %s (%i) (%s)\n", topIn.TruncResNameOnumId(*rnum).c_str(), topIn.Res(*rnum).OriginalResNum(), topIn.Res(*rnum).chainID());
    int nd1idx = -1;
    int ne2idx = -1;
    Residue const& hisRes = topIn.Res( *rnum );
    for (int at = hisRes.FirstAtom(); at < hisRes.LastAtom(); ++at)
    {
      if ( (topIn[at].Name() == ND1 ) )
        nd1idx = at;
      else if ( (topIn[at].Name() == NE2 ) )
        ne2idx = at;
    }
    if (nd1idx == -1) {
      mprintf("Warning: Atom %s not found for %s; skipping residue.\n", *ND1, topIn.TruncResNameOnumId(*rnum).c_str());
      continue;
    }
    if (ne2idx == -1) {
      mprintf("Warning: Atom %s not found for %s; skipping residue,\n", *NE2, topIn.TruncResNameOnumId(*rnum).c_str());
      continue;
    }
    if (StructureDebugLevel() > 1)
      mprintf("DEBUG: %s nd1idx= %i ne2idx= %i\n",
              topIn.TruncResNameOnumId( *rnum ).c_str(), nd1idx+1, ne2idx+1);
    // Check for H bonded to nd1/ne2
    int nd1h = 0;
    for (Atom::bond_iterator bat = topIn[nd1idx].bondbegin();
                             bat != topIn[nd1idx].bondend();
                           ++bat)
      if ( topIn[*bat].Element() == Atom::HYDROGEN)
        ++nd1h;
    if (nd1h > 1) {
      mprinterr("Error: More than 1 hydrogen bonded to %s\n",
                topIn.ResNameNumAtomNameNum(nd1idx).c_str());
      return 1;
    }
    int ne2h = 0;
    for (Atom::bond_iterator bat = topIn[ne2idx].bondbegin();
                             bat != topIn[ne2idx].bondend();
                           ++bat)
      if ( topIn[*bat].Element() == Atom::HYDROGEN)
        ++ne2h;
    if (ne2h > 1) {
      mprinterr("Error: More than 1 hydrogen bonded to %s\n",
                topIn.ResNameNumAtomNameNum(ne2idx).c_str());
      return 1;
    }
    if (nd1h > 0 && ne2h > 0) {
      mprintf("\t\t%s => %s\n", topIn.TruncResNameOnumId(*rnum).c_str(), *HipName);
      ChangeResName( topIn.SetRes(*rnum), HipName );
      nchanged++;
    } else if (nd1h > 0) {
      mprintf("\t\t%s => %s\n", topIn.TruncResNameOnumId(*rnum).c_str(), *HidName);
      ChangeResName( topIn.SetRes(*rnum), HidName );
      nchanged++;
    } else if (ne2h > 0) {
      mprintf("\t\t%s => %s\n", topIn.TruncResNameOnumId(*rnum).c_str(), *HieName);
      ChangeResName( topIn.SetRes(*rnum), HieName );
      nchanged++;
    }
    //else {
    //  // Default to epsilon
    //  mprintf("\tUsing default name '%s' for %s\n", *HieName, topIn.TruncResNameOnumId(*rnum).c_str());
    //  HisResNames.push_back( HieName );
    //}
  }
  if (nchanged == 0) 
    mprintf("\tNo histidine names were changed.\n");
  else
    mprintf("\t%u histidine names were changed.\n", nchanged);
  return 0;
}


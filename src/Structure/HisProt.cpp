#include "HisProt.h"
#include "StructureRoutines.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../NameType.h"
#include "../Topology.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
HisProt::HisProt() {}

/// Warn if name length greater than NameType max
static inline void warn_length(std::string const& nameIn) {
  if (nameIn.size() >= NameType::max())
    mprintf("Warning: Name '%s' is too large and will be truncated.\n", nameIn.c_str());
}

const char* HisProt::keywords_ =
  "\t[{nohisdetect |\n"
  "\t  [nd1 <nd1>] [ne2 <ne2] [hisname <his>] [hiename <hie>]\n"
  "\t  [hidname <hid>] [hipname <hip]} [defaulthis <default>]]\n";

/** Initialize from args. Get atom/residue names. */
int HisProt::InitHisProt(ArgList& argIn, int debugIn) {
  nd1name_ = argIn.GetStringKey("nd1", "ND1");
  ne2name_ = argIn.GetStringKey("ne2", "NE2");
  hisname_ = argIn.GetStringKey("hisname", "HIS");
  hiename_ = argIn.GetStringKey("hiename", "HIE");
  hidname_ = argIn.GetStringKey("hidname", "HID");
  hipname_ = argIn.GetStringKey("hipname", "HIP");
  default_ = argIn.GetStringKey("defaulthis");
  warn_length(nd1name_);
  warn_length(ne2name_);
  warn_length(hisname_);
  warn_length(hiename_);
  warn_length(hidname_);
  warn_length(hipname_);
  if (!default_.empty())
    warn_length(default_);
  return 0;
}

/** Print info to stdout. */
void HisProt::HisProtInfo() const {
  mprintf("\tHistidine protonation detection:\n");
  mprintf("\t\tND1 atom name                   : %s\n", nd1name_.c_str());
  mprintf("\t\tNE2 atom name                   : %s\n", ne2name_.c_str());
  mprintf("\t\tHistidine original residue name : %s\n", hisname_.c_str());
  mprintf("\t\tEpsilon-protonated residue name : %s\n", hiename_.c_str());
  mprintf("\t\tDelta-protonated residue name   : %s\n", hidname_.c_str());
  mprintf("\t\tDoubly-protonated residue name  : %s\n", hipname_.c_str());
  if (!default_.empty())
     mprintf("\t\tDefault residue name            : %s\n", default_.c_str());
}

/** Determine histidine protonation from hydrogens. */
int HisProt::DetermineHisProt(Topology& topIn) const {
  return determineHisProt(topIn, nd1name_, ne2name_, hisname_,
                          hiename_, hiename_, hipname_, default_);
}

/** Try to determine histidine protonation from existing hydrogens.
  * Change residue names as appropriate.
  * \param topIn Input topology.
  * \param ND1 ND1 atom name.
  * \param NE2 NE2 atom name.
  * \param HisName PDB histidine name.
  * \param HieName Name for epsilon-protonated His.
  * \param HidName Name for delta-protonated His.
  * \param HipName Name for doubly-protonated His.
  * \param defaultName optional residue name for default if nothing else detected
  */
int HisProt::determineHisProt( Topology& topIn,
                                           NameType const& ND1,
                                           NameType const& NE2,
                                           NameType const& HisName,
                                           NameType const& HieName,
                                           NameType const& HidName,
                                           NameType const& HipName,
                                           std::string const& defaultName )
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
    } else if (!defaultName.empty()) {
      // Default
      NameType defName(defaultName);
      mprintf("\tUsing default name '%s' for %s\n", *defName, topIn.TruncResNameOnumId(*rnum).c_str());
      ChangeResName( topIn.SetRes(*rnum), defName );
      nchanged++;
    }
  }
  if (nchanged == 0) 
    mprintf("\tNo histidine names were changed.\n");
  else
    mprintf("\t%u histidine names were changed.\n", nchanged);
  return 0;
}


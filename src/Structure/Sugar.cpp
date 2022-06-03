#include "Sugar.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"

using namespace Cpptraj::Structure;

/** CONSTRUCTOR - Incomplete setup; set anomeric atom as residue first atom
  *               so that ResNum() works.
  */
Sugar::Sugar(StatType status, int firstat) :
  stat_(status),
  ring_oxygen_atom_(-1),
  anomeric_atom_(firstat),
  ano_ref_atom_(-1),
  highest_stereocenter_(-1),
  ringType_(SugarToken::UNKNOWN_RING),
  isMissingAtoms_(false)
{}

/** CONSTRUCTOR - Partial setup; ring O, anomeric atom, ring and chain atoms. */
Sugar::Sugar(StatType status, int roa, int aa,
                                  Iarray const& RA, Iarray const& CA) :
  stat_(status),
  ring_oxygen_atom_(roa),
  anomeric_atom_(aa),
  ano_ref_atom_(-1),
  highest_stereocenter_(-1),
  ringType_(SugarToken::UNKNOWN_RING),
  isMissingAtoms_(true),
  ring_atoms_(RA),
  chain_atoms_(CA)
{
  if (ring_oxygen_atom_ != -1) {
    if (RA.size() == 5)
      ringType_ = SugarToken::PYRANOSE;
    else if (RA.size() == 4)
      ringType_ = SugarToken::FURANOSE;
  }
}

/** CONSTRUCTOR - Set ring atom indices and ring type. */
Sugar::Sugar(int roa, int aa, int ara, int hs,
                                  Iarray const& RA, Iarray const& CA, bool isMissing) :
  stat_(SETUP_OK),
  ring_oxygen_atom_(roa),
  anomeric_atom_(aa),
  ano_ref_atom_(ara),
  highest_stereocenter_(hs),
  ringType_(SugarToken::UNKNOWN_RING),
  isMissingAtoms_(isMissing),
  ring_atoms_(RA),
  chain_atoms_(CA)
{
  if (ring_oxygen_atom_ != -1) {
    if (RA.size() == 5)
      ringType_ = SugarToken::PYRANOSE;
    else if (RA.size() == 4)
      ringType_ = SugarToken::FURANOSE;
  }
}

/** \return Residue number based on the anomeric carbon atom index. */
int Sugar::ResNum(Topology const& topIn) const {
  return topIn[anomeric_atom_].ResNum();
}

/** Correspond to Sugar::StatType */
const char* Sugar::StatTypeStr_[] = {
  "OK", "Missing ring oxygen", "Multiple ring oxygens", "Missing chain atoms",
  "Missing anomeric reference", "Missing configurational carbon" };

/** Print info about the sugar to STDOUT. */
void Sugar::PrintInfo(Topology const& topIn) const {
  if (NotSet()) {
    mprintf("\t%s : Not Set. %s\n", topIn.TruncResNameOnumId(ResNum(topIn)).c_str(),
            StatTypeStr_[stat_]);
  } else {
    mprintf("\t%s : %s\n", topIn.TruncResNameOnumId(ResNum(topIn)).c_str(),
            StatTypeStr_[stat_]);
  }
  //mprintf("\t%s :\n", topIn.TruncResNameOnumId(ResNum(topIn)).c_str());
  if (isMissingAtoms_)
    mprintf("\t\tIs missing atoms.\n");
  if (ring_oxygen_atom_ != -1)
    mprintf("\t\tRing O           : %s\n", topIn.TruncAtomNameNum(ring_oxygen_atom_).c_str());
  if (anomeric_atom_ != -1)
    mprintf("\t\tAnomeric C       : %s\n", topIn.TruncAtomNameNum(anomeric_atom_).c_str());
  if (ano_ref_atom_ != -1)
    mprintf("\t\tAnomeric ref. C  : %s\n", topIn.TruncAtomNameNum(ano_ref_atom_).c_str());
  if (highest_stereocenter_ != -1)
    mprintf("\t\tConfig. C        : %s\n", topIn.TruncAtomNameNum(highest_stereocenter_).c_str());
  mprintf("\t\tNum ring atoms   : %u\n", NumRingAtoms());
  static const char* RingPrefixStr[] = { "Hexo", "Pento", "Other" };
  int rpsidx = 2;
  if (chain_atoms_.size() == 6)
    rpsidx = 0;
  else if (chain_atoms_.size() == 5)
    rpsidx = 1;
  static const char* RingTypeStr[] = { "Pyranose", "Furanose", "Unknown" };
  mprintf("\t\tRing Type        : %s%s\n", RingPrefixStr[rpsidx], RingTypeStr[ringType_]);
  mprintf("\t\tNon-O Ring atoms :");
  for (Iarray::const_iterator it = ring_atoms_.begin(); it != ring_atoms_.end(); ++it)
    mprintf(" %s", topIn.TruncAtomNameNum(*it).c_str());
  mprintf("\n\t\tChain atoms      :");
  for (Iarray::const_iterator it = chain_atoms_.begin(); it != chain_atoms_.end(); ++it)
    mprintf(" %s", topIn.TruncAtomNameNum(*it).c_str());
  mprintf("\n");
}

/** \return Number of ring atoms (including oxygen) */
unsigned int Sugar::NumRingAtoms() const {
  if (NotSet()) return 0;
  return ring_atoms_.size() + 1;
}

/// \return what oldidx should be according to atomMap
inline static int find_new_idx(int oldidx, std::vector<int> const& atomMap, int at0, int at1) {
  for (int newidx = at0; newidx != at1; newidx++)
    if (atomMap[newidx] == oldidx)
      return newidx;
  return -1;
}

/** Remap internal indices according to given map. */
void Sugar::RemapIndices(Iarray const& atomMap, int at0, int at1) {
  // Always try the anomeric atom
  anomeric_atom_ = find_new_idx(anomeric_atom_, atomMap, at0, at1);
  //mprintf("DEBUG: Ring O old = %i ring O new %i\n", ring_oxygen_atom_+1, atomMap[ring_oxygen_atom_]+1);
  if (ring_oxygen_atom_ != -1)
    ring_oxygen_atom_     = find_new_idx(ring_oxygen_atom_, atomMap, at0, at1);
  if (ano_ref_atom_ != -1)
    ano_ref_atom_         = find_new_idx(ano_ref_atom_, atomMap, at0, at1);
  if (highest_stereocenter_ != -1)
    highest_stereocenter_ = find_new_idx(highest_stereocenter_, atomMap, at0, at1);
  for (Iarray::iterator it = ring_atoms_.begin(); it != ring_atoms_.end(); ++it)
    *it = find_new_idx(*it, atomMap, at0, at1);
  for (Iarray::iterator it = chain_atoms_.begin(); it != chain_atoms_.end(); ++it)
    *it = find_new_idx(*it, atomMap, at0, at1);
}



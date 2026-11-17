#include "ModXNA_Info.h"
#include "ArgList.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // SplitString

/** CONSTRUCTOR */
ModXNA_Info::ModXNA_Info() {}

/** Array of modxna mask types for the header */
const char* ModXNA_Info::MaskTypeStr_[] = {
  "HEAD01", "HEAD01STRIP", "TAIL01", "TAIL01STRIP", "HEAD03", "HEAD03STRIP"
};

/** Given a ModXNA metadata string, return an array of just the
  * atom names.
  */
ModXNA_Info::Sarray ModXNA_Info::atomNameList(std::string const& mheader, int& errval)
{
  Sarray ListOut;
  if (mheader.empty()) return ListOut;
  // Use SplitString instead of ArgList to handle atom names
  // with leading quotes.
  Sarray anames = SplitString(mheader, "@,");
  bool has_none = false;
  for (unsigned int iarg = 0; iarg < anames.size(); iarg++) {
    //mprintf("DEBUG: Atom name %u : %s\n", iarg, anames[iarg].c_str());
    // NONE is a special case of blank atom name. TODO check everywhere
    if ( anames[iarg] == "NONE" )
      has_none = true;
    else
      ListOut.push_back( anames[iarg] );
  }
  if (ListOut.empty() && !has_none) {
    mprinterr("Error: ModXNA header string '%s' has no valid atom names.\n", mheader.c_str());
    errval++;
  }
  return ListOut;
}

/** Parse ModXNA info from given string. */
ModXNA_Info::StatType ModXNA_Info::ParseModxnaStr(std::string const& line) {
  if (line.empty()) return NOT_MODXNA;
  ArgList args(line, ":");
  //args.PrintDebug(); // DEBUG

  if (args.Nargs() < 1)
    return NOT_MODXNA;
  else if (args[0].compare(0,7,"modXNA-") != 0)
    return NOT_MODXNA;
  modxnaType_ = args[0];
  args.MarkArg(0);

  if (args.Nargs() < 2) {
    mprinterr("Error: modXNA string is too short.\n"
              "Error: '%s'\n", line.c_str());
    return ERR;
  }
  fragmentName_ = args[1];
  args.MarkArg(1);

  int errval = 0;
  head_ = atomNameList(args.GetStringKey("HEAD01"), errval);
  headStrip_ = atomNameList(args.GetStringKey("HEAD01STRIP"), errval);
  tail_ = atomNameList(args.GetStringKey("TAIL01"), errval);
  tailStrip_ = atomNameList(args.GetStringKey("TAIL01STRIP"), errval);
  anchor_ = atomNameList(args.GetStringKey("HEAD03"), errval);
  anchorStrip_ = atomNameList(args.GetStringKey("HEAD03STRIP"), errval);

  if (errval > 0) return ERR;
  return OK;
}

/** Generate mask string from atom names. */
std::string ModXNA_Info::atomMaskStr(Sarray const& anames, std::string const& noMask)
{
  if (anames.empty()) return noMask;
  std::string maskstr("@" + anames[0]);
  for (unsigned int idx = 1; idx < anames.size(); idx++)
    maskstr.append("," + anames[idx]);
  return maskstr;
}

/** Generate mask string from atom names. */
std::string ModXNA_Info::atomMaskStr(Sarray const& anames) {
  return atomMaskStr(anames, "");
}

/** Print name array to stdout. */
void ModXNA_Info::printmxna(const char* desc, Sarray const& anames) {
  if (!anames.empty()) {
    std::string str = atomMaskStr(anames);
    if (!str.empty()) mprintf("%s: %s\n", desc, str.c_str());
  }
}

/** Print info to stdout */
void ModXNA_Info::PrintModxna() const {
  if (!HasModxna()) return;
  mprintf("\tType        : %s\n", modxnaType_.c_str());
  mprintf("\tFragment    : %s\n", fragmentName_.c_str());
  printmxna("\tHead        ", head_);
  printmxna("\tHeadStrip   ", headStrip_);
  printmxna("\tTail        ", tail_);
  printmxna("\tTailStrip   ", tailStrip_);
  printmxna("\tAnchor      ", anchor_);
  printmxna("\tAnchorStrip ", anchorStrip_);
}

/** \return specified array of names */
ModXNA_Info::Sarray const& ModXNA_Info::AtomNames(MaskType typeIn) const {
  switch (typeIn) {
    case HEAD        : return head_;
    case HEADSTRIP   : return headStrip_;
    case TAIL        : return tail_;
    case TAILSTRIP   : return tailStrip_;
    case ANCHOR      : return anchor_;
    case ANCHORSTRIP : return anchorStrip_;
  }
  return head_;
}

/** \return Atom mask string for the specified mask type */
std::string ModXNA_Info::AtomMaskString(MaskType typeIn) const {
  switch (typeIn) {
    case HEAD        : return atomMaskStr(head_);
    case HEADSTRIP   : return atomMaskStr(headStrip_);
    case TAIL        : return atomMaskStr(tail_);
    case TAILSTRIP   : return atomMaskStr(tailStrip_);
    case ANCHOR      : return atomMaskStr(anchor_);
    case ANCHORSTRIP : return atomMaskStr(anchorStrip_);
  }
  return std::string("");
}

/**  \return ModXNA mask type containing the given atom name, -1 if not present.
  */
int ModXNA_Info::FindMaskWithAtom(std::string const& nameIn) const {
  for (int idx = 0; idx <= (int)ANCHORSTRIP; idx++) {
    Sarray const& anames = AtomNames( (MaskType)idx );
    for (Sarray::const_iterator it = anames.begin(); it != anames.end(); ++it) {
      if (*it == nameIn) return idx;
    }
  }
  return -1;
}

static inline void rename_atoms(std::vector<std::string>& anames,
                                std::string const& oldname, std::string const& newname)
{
  for (std::vector<std::string>::iterator it = anames.begin(); it != anames.end(); ++it) {
    if (*it == oldname)
      *it = newname;
  }
}

/** Rename any atoms with given name to new name. */
void ModXNA_Info::RenameAtomInMasks(std::string const& oldname, std::string const& newname) {
  rename_atoms( head_, oldname, newname );
  rename_atoms( headStrip_, oldname, newname );
  rename_atoms( tail_, oldname, newname );
  rename_atoms( tailStrip_, oldname, newname );
  rename_atoms( anchor_, oldname, newname );
  rename_atoms( anchorStrip_, oldname, newname );
}

/** Generate ModXNA metadata string */
std::string ModXNA_Info::GenMetadataString() const {
  if (fragmentName_.empty()) return std::string("");
  std::string out(modxnaType_ + ":" + fragmentName_);
  // Base
  out.append(":HEAD01:" + atomMaskStr(head_, "NONE"));
  out.append(":HEAD01STRIP:" + atomMaskStr(headStrip_, "NONE"));
  if (!tail_.empty()) {
    // BB
    out.append(":TAIL01:" + atomMaskStr(tail_, "NONE"));
    out.append(":TAIL01STRIP:" + atomMaskStr(tailStrip_, "NONE"));
    if (!anchor_.empty()) {
      // Sugar
      out.append(":HEAD03:" + atomMaskStr(anchor_, "NONE"));
      out.append(":HEAD03STRIP:" + atomMaskStr(anchorStrip_, "NONE"));
    }
  }
  return out;
}



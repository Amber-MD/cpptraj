#include "ModXNA_Info.h"
#include "ArgList.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
ModXNA_Info::ModXNA_Info() {}

/** Parse ModXNA info from given string. */
ModXNA_Info::StatType ModXNA_Info::ParseModxnaStr(std::string const& line) {
  if (line.empty()) return NOT_MODXNA;
  ArgList args(line, ":");

  if (args.GetStringNext() != "modXNA-fragment") return NOT_MODXNA;

  if (args.Nargs() < 2) {
    mprinterr("Error: modXNA string is too short.\n"
              "Error: '%s'\n", line.c_str());
    return ERR;
  }
  fragmentName_ = args[1];

  if (args.Nargs() >= 4) head_ = args[3];
  if (args.Nargs() >= 6) headStrip_ = args[5];
  if (args.Nargs() >= 8) tail_ = args[7];
  if (args.Nargs() >= 10) tailStrip_ = args[9];
  if (args.Nargs() >= 12) anchor_ = args[11];
  if (args.Nargs() >= 14) anchorStrip_ = args[13];

  return OK;
}

static inline void printmxna(const char* desc, std::string const& str) {
  if (!str.empty()) mprintf("%s: %s\n", desc, str.c_str());
}

/** Print info to stdout */
void ModXNA_Info::PrintModxna() const {
  if (!HasModxna()) return;
  mprintf("\tFragment: %s\n", fragmentName_.c_str());
  printmxna("\tHead        ", head_);
  printmxna("\tHeadStrip   ", headStrip_);
  printmxna("\tTail        ", tail_);
  printmxna("\tTailStrip   ", tailStrip_);
  printmxna("\tAnchor      ", anchor_);
  printmxna("\tAnchorStrip ", anchorStrip_);
}

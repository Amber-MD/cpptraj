#include "EnsembleOut.h"
#include "CpptrajStdio.h"

Range EnsembleOut::MembersToWrite(std::string const& onlyMembers, int ensembleSize) {
  Range members;
  int err;
  // Empty String indicates write All members
  if (onlyMembers.empty())
    err = members.SetRange(0, ensembleSize);
  else
    err = members.SetRange(onlyMembers);
  if (err != 0 || members.Empty()) {
    mprinterr("Error: onlymembers: Invalid range (%s)\n", onlyMembers.c_str());
    return Range();
  }
  return members;
}

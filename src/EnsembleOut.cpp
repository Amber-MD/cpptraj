#include "EnsembleOut.h"
#include "CpptrajStdio.h"

int EnsembleOut::SetMembersToWrite(std::string const& onlyMembers, int ensembleSize) {
  members_to_write_.Clear();
  int err;
  // Empty String indicates write All members
  if (onlyMembers.empty())
    err = members_to_write_.SetRange(0, ensembleSize);
  else
    err = members_to_write_.SetRange(onlyMembers);
  if (err != 0 || members_to_write_.Empty()) {
    mprinterr("Error: onlymembers: Invalid range (%s)\n", onlyMembers.c_str());
    return 1;
  }
  return 0;
}

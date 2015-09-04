#include "InputTrajCommon.h"
#include "CpptrajStdio.h"

int InputTrajCommon::SetNameAndParm(std::string const& fname, Topology* top) {
  if (top == 0) {
    mprinterr("Internal Error: Trajin::SetNameAndParm(): Topology is null.\n");
    return 1;
  }
  trajParm_ = top;
  if (fname.empty()) {
    mprinterr("Internal Error: Trajin::SetNameAndParm(): File name is empty.\n");
    return 1;
  }
  if (trajName_.SetFileNameWithExpansion( fname )) return 1;
  return 0;
}

#include "InputTrajCommon.h"
#include "CpptrajStdio.h"

int InputTrajCommon::SetNameAndParm(FileName const& fname, Topology* top) {
  if (top == 0) {
    mprinterr("Internal Error: Trajin::SetNameAndParm(): Topology is null.\n");
    return 1;
  }
  trajParm_ = top;
  if (fname.empty()) {
    mprinterr("Internal Error: Trajin::SetNameAndParm(): File name is empty.\n");
    return 1;
  }
  trajName_ = fname;
  if (!File::Exists( trajName_ )) {
    File::ErrorMsg( trajName_.full() );
    return 1;
  }
  return 0;
}

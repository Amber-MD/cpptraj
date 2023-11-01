#include "AssociatedData_Connect.h"
#include "CpptrajStdio.h"

const char* AssociatedData_Connect::HelpText = "";

int AssociatedData_Connect::ProcessAdataArgs(ArgList&) {
  return 0;
}

void AssociatedData_Connect::Ainfo() const {
  if (!connect_.empty()) {
    mprintf(" (Connect Atoms:");
    for (Iarray::const_iterator it = connect_.begin(); it != connect_.end(); ++it)
      mprintf(" %i", *it + 1);
    mprintf(")");
  }
  return;
}

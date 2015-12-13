#include <cstdlib>
#include "Exec_System.h"
#include "CpptrajStdio.h"

void Exec_System::Help() const { mprintf("  Call command from system.\n"); }

Exec::RetType Exec_System::Execute(CpptrajState& State, ArgList& argIn) {
  int err = system( argIn.ArgLine() );
  if (err != 0) mprintf("Warning: '%s' returned %i\n", argIn.Command(), err);
  return CpptrajState::OK;
}

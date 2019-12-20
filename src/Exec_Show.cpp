#include "Exec_Show.h"
#include "CpptrajStdio.h"

void Exec_Show::Help() const {
  mprintf("  Show all current script variables and their values.\n");
}

Exec::RetType Exec_Show::Execute(CpptrajState& State, ArgList& argIn)
{
  State.DSL().ListStringVar();
  return CpptrajState::OK;
}

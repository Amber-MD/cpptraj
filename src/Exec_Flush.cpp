#include "Exec_Flush.h"
#include "CpptrajStdio.h"

// Exec_Flush::Help()
void Exec_Flush::Help() const
{
  mprintf("  Write any pending data files; close all other files.\n");
}

// Exec_Flush::Execute()
Exec::RetType Exec_Flush::Execute(CpptrajState& State, ArgList& argIn)
{
  State.DFL().Flush();
  return CpptrajState::OK;
}

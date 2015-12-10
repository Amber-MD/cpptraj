#include "Exec_Commands.h"
#include "CpptrajStdio.h"

void Exec_Run::Help() const {
  mprintf("  Process all trajectories currently in input trajectory list.\n"
          "  All actions in action list will be run on each frame.\n"
          "  If not processing ensemble input, all analyses in analysis\n"
          "  list will be run after trajectory processing.\n");
}

void Exec_NoExitOnError::Help() const {
  mprintf("  Do not exit when errors are encountered. This is the default\n"
          "  in interactive mode.\n");
}

Exec::RetType Exec_NoExitOnError::Execute(CpptrajState& State, ArgList&)
{
  State.SetNoExitOnError();
  mprintf("\tAttempting to ignore errors if possible.\n");
  return CpptrajState::OK;
}

void Exec_NoProgress::Help() const {
  mprintf("  Do not print progress while reading in trajectories.\n");
}

Exec::RetType Exec_NoProgress::Execute(CpptrajState& State, ArgList&)
{
  State.SetNoProgress();
  mprintf("\tProgress bar will not be used during Run.\n");
  return CpptrajState::OK;
}

void Exec_Quit::Help() const { mprintf("  Exit CPPTRAJ\n"); }

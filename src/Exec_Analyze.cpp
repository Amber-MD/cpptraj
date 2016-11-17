#include "Exec_Analyze.h"
#include "CpptrajStdio.h"
#include "Command.h"

// Exec_Analyze::Help()
void Exec_Analyze::Help() const {
  mprintf("\t[<analysis command>]\n"
          "  Add specified analysis command to analysis queue.\n"
          "  This is retained for backwards compatibility only. Users are\n"
          "  encouraged to use analysis commands themselves.\n");
}

// Exec_Analyze::Execute()
Exec::RetType Exec_Analyze::Execute(CpptrajState& State, ArgList& argIn) {
  // Remove 'analyze'
  ArgList arg = argIn;
  arg.RemoveFirstArg();
  if (arg.empty()) {
    mprinterr("Error: No analysis command specified.\n");
    return CpptrajState::ERR;
  }
  Cmd const& cmd = Command::SearchTokenType(DispatchObject::ANALYSIS, arg.Command());
  if (cmd.Empty()) {
    mprinterr("Error: Analysis command '%s' not found.\n", arg.Command());
    return CpptrajState::ERR;
  }
  return State.AddToAnalysisQueue( (Analysis*)cmd.Alloc(), arg );
}

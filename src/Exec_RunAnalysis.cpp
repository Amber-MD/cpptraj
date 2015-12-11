#include "Exec_RunAnalysis.h"
#include "CpptrajStdio.h"
#include "Command.h"
#include "Timer.h"

void Exec_RunAnalysis::Help() const {
  mprintf("\t[<analysis> [<analysis args>]]\n"
          "  If specified alone, run all analyses in the analysis list.\n"
          "  Otherwise run the specified analysis immediately.\n");
}

Exec::RetType Exec_RunAnalysis::Execute(CpptrajState& State, ArgList& argIn) {
  // If only 1 arg (the command) run all analyses in list
  if (argIn.Nargs() == 1) {
    int eval = State.RunAnalyses();
    State.MasterDataFileWrite();
    if (eval == 0)
      return CpptrajState::OK;
    else
      return CpptrajState::ERR;
  }
  // Run specified analysis
  // FIXME: Use RemoveFirstArg
  ArgList analyzeargs = argIn.RemainingArgs();
  analyzeargs.MarkArg(0);
  Cmd const& cmd = Command::SearchTokenType( DispatchObject::ANALYSIS, analyzeargs.Command() );
  if ( cmd.Empty() ) return CpptrajState::ERR;
  Analysis* ana = (Analysis*)cmd.Alloc();
  if (ana == 0) return CpptrajState::ERR;
  Timer total_time;
  total_time.Start();
  CpptrajState::RetType err = CpptrajState::ERR;
  AnalysisSetup setup(*State.DSL(), *State.DFL());
  if ( ana->Setup( analyzeargs, setup, State.Debug() ) == Analysis::OK )
  {
    analyzeargs.CheckForMoreArgs();
    if (ana->Analyze() != Analysis::ERR) {
      err = CpptrajState::OK;
      State.MasterDataFileWrite();
    }
  }
  delete ana;
  total_time.Stop();
  mprintf("TIME: Total analysis execution time: %.4f seconds.\n", total_time.Total());
  return err;
}

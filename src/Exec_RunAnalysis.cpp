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
  int err = 0;
# ifdef MPI
  // Only master performs analyses currently.
  if (Parallel::TrajComm().Size() > 1)
    mprintf("Warning: Analysis does not currently use multiple MPI threads.\n");
  if (Parallel::TrajComm().Master())
# endif
    err = DoRunAnalysis(State, argIn);
# ifdef MPI
  if (Parallel::World().CheckError( err )) err = 1;
# endif
  if (err != 0) return CpptrajState::ERR;
  return CpptrajState::OK;
}

int Exec_RunAnalysis::DoRunAnalysis(CpptrajState& State, ArgList& argIn) const {
  // FIXME: Use RemoveFirstArg
  ArgList analyzeargs = argIn.RemainingArgs();
  analyzeargs.MarkArg(0);
  Cmd const& cmd = Command::SearchTokenType( DispatchObject::ANALYSIS, analyzeargs.Command() );
  if ( cmd.Empty() ) return 1;
  Analysis* ana = (Analysis*)cmd.Alloc();
  if (ana == 0) return 1;
  Timer total_time;
  total_time.Start();
  CpptrajState::RetType err = CpptrajState::ERR;
  AnalysisSetup setup(State.DSL(), State.DFL());
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
  if (err != CpptrajState::OK) return 1;
  return 0;
}

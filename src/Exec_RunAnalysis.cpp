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
  //rprintf("DEBUG: Entering runanalysis. Nargs is %i\n", argIn.Nargs());
  int err = 0;
  if (argIn.Nargs() == 1) {
    // If only 1 arg (the command) run all analyses in list
    err = State.RunAnalyses();
  } else {
    // Run specified analysis
    err = DoRunAnalysis(State, argIn);
  }
  State.MasterDataFileWrite();
# ifdef MPI
  if (Parallel::World().CheckError( err )) err = 1;
# endif
  if (err != 0) return CpptrajState::ERR;
  return CpptrajState::OK;
}

/** Set up and perform the analysis. */
int Exec_RunAnalysis::DoRunAnalysis(CpptrajState& State, ArgList& argIn) const {
  //rprintf("DEBUG: Entering DoRunAnalysis.\n");
  // FIXME: Use RemoveFirstArg
  ArgList analyzeargs = argIn.RemainingArgs();
  analyzeargs.MarkArg(0);
  Cmd const& cmd = Command::SearchTokenType( DispatchObject::ANALYSIS, analyzeargs.Command() );
  //rprintf("DEBUG: cmd is empty= %i\n", (int)cmd.Empty());
  if ( cmd.Empty() ) return 1;
  Analysis* ana = (Analysis*)cmd.Alloc();
  if (ana == 0) return 1;
  //rprintf("DEBUG: Starting analysis.\n");

  Timer total_time;
  total_time.Start();

  bool run_analysis = true;
# ifdef MPI
  Parallel::Comm const& analyzeComm = Parallel::World();
  AnalysisSetup setup(State.DSL(), State.DFL(), analyzeComm);
  if (!ana->IsParallel()) {
    mprintf("Warning: Analysis '%s' does not currently use multiple MPI processes.\n", analyzeargs.Command());
    run_analysis = analyzeComm.Master();
  }
# else /* MPI */
  AnalysisSetup setup(State.DSL(), State.DFL());
# endif /* MPI */
  Analysis::RetType ret = Analysis::OK;
  //rprintf("DEBUG: Run analysis= %i\n",(int)run_analysis);
  if (run_analysis) {
    ret = ana->Setup( analyzeargs, setup, State.Debug() );
    if (analyzeargs.CheckForMoreArgs()) {
      ret = Analysis::ERR;
    } else if (ret == Analysis::OK) {
      ret = ana->Analyze();
    }
  }
  int stat;
  if (ret == Analysis::ERR)
    stat = 1;
  else
    stat = 0;

  delete ana;

  total_time.Stop();
  mprintf("TIME: Total analysis execution time: %.4f seconds.\n", total_time.Total());
  return stat;
}

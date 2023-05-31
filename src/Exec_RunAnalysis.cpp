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
  int err = DoRunAnalysis(State, argIn);
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
  CpptrajState::RetType stat = CpptrajState::ERR;
# ifdef MPI
  AnalysisSetup setup(State.DSL(), State.DFL(), Parallel::World());
# else
  AnalysisSetup setup(State.DSL(), State.DFL());
# endif
  if ( ana->Setup( analyzeargs, setup, State.Debug() ) == Analysis::OK )
  {
    analyzeargs.CheckForMoreArgs();
    Analysis::RetType ret;
#   ifdef MPI
    if (ana->IsParallel()) {
      ret = ana->Analyze();
    } else {
      mprintf("Warning: Analysis '%s' does not currently use multiple MPI processes.\n", analyzeargs.Command());
      if ( Parallel::World().Master() )
        ret = ana->Analyze();
      Parallel::World().MasterBcast( &ret, 1, MPI_INT );
    }
    int err;
    if (ret == Analysis::ERR) {
      rprinterr("Error: In parallel, analysis '%s' failed.\n", analyzeargs.Command());
      err = 1;
    } else
      err = 0;
    if (Parallel::World().CheckError( err ))
      ret = Analysis::ERR;
#   else /* MPI */
    ret = ana->Analyze();
#   endif /* MPI */
    if (ret != Analysis::ERR) {
      //rprintf("DEBUG: Analysis success!\n");
      stat = CpptrajState::OK;
      State.MasterDataFileWrite();
    }
  }
  delete ana;
  total_time.Stop();
  mprintf("TIME: Total analysis execution time: %.4f seconds.\n", total_time.Total());
  if (stat != CpptrajState::OK) return 1;
  return 0;
}

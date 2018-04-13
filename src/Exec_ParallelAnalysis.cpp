#include "Exec_ParallelAnalysis.h"
#include "CpptrajStdio.h"
#ifdef MPI
#include <cstdio> // DEBUG
// Exec_ParallelAnalysis::Help()
void Exec_ParallelAnalysis::Help() const
{

}

// Exec_ParallelAnalysis::Execute()
Exec::RetType Exec_ParallelAnalysis::Execute(CpptrajState& State, ArgList& argIn)
{
  // DEBUG - Have each thread report what analyses it knows about and what
  // data sets it has.
/*
  for (int rank = 0; rank < Parallel::World().Size(); rank++) {
    if (rank == Parallel::World().Rank()) {
      printf("Rank %i, %u analyses, %zu data sets:\n", rank, State.Analyses().size(),
             State.DSL().size());
      for (DataSetList::const_iterator it = State.DSL().begin(); it != State.DSL().end(); ++it)
        printf("\t'%s' (%zu)\n", (*it)->Meta().PrintName().c_str(), (*it)->Size());
    }
    Parallel::World().Barrier();
  }
  Parallel::World().Barrier();
*/
  // Naively divide up all analyses among threads.
  int my_start, my_stop;
  int nelts = Parallel::World().DivideAmongThreads( my_start, my_stop, State.Analyses().size() );
  rprintf("Dividing %zu analyses among %i threads: %i to %i (%i)\n",
          State.Analyses().size(), Parallel::World().Size(), my_start, my_stop, nelts);

  // Each thread runs the analyses they are responsible for.
  int nerr = 0;
  for (int na = my_start; na != my_stop; na++) {
    // TODO check setup status
    if (State.Analyses().Ana(na).Analyze() != Analysis::OK) {
      rprinterr("Error: Analysis failed: '%s'\n", State.Analyses().Args(na).ArgLine());
      nerr++;
    }
  }
  // This error check serves as a barrier
  if (Parallel::World().CheckError( nerr )) return CpptrajState::ERR;
  State.Analyses().Clear();
  return CpptrajState::OK;
}
#else
// Exec_ParallelAnalysis::Help()
void Exec_ParallelAnalysis::Help() const
{
  mprintf("  This command is only available in MPI-enabled builds.\n");
}

// Exec_ParallelAnalysis::Execute()
Exec::RetType Exec_ParallelAnalysis::Execute(CpptrajState& State, ArgList& argIn)
{
  mprinterr("Error: This command is only available in MPI-enabled builds.\n");
  return CpptrajState::ERR;
}
#endif

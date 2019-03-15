#include <vector>
#include "Exec_ParallelAnalysis.h"
#include "CpptrajStdio.h"
#ifdef MPI
#include "Timer.h"

// Exec_ParallelAnalysis::Help()
void Exec_ParallelAnalysis::Help() const
{
  mprintf("\t[sync]\n"
          "  Divide all currently set up Analyses among available MPI processes and run.\n"
          "  If 'sync' is specifed, sync all results back to master process; this may be\n"
          " required if performing subsequent analyses.\n");
}

// Exec_ParallelAnalysis::Execute()
Exec::RetType Exec_ParallelAnalysis::Execute(CpptrajState& State, ArgList& argIn)
{
  Timer t_total;
  t_total.Start();
  Timer t_sync;
  bool syncToMaster = argIn.hasKey("sync");
  std::vector<unsigned int> setSizesBefore;
  if (syncToMaster) {
    t_sync.Start();
    setSizesBefore.reserve( State.DSL().size() );
    for (DataSetList::const_iterator it = State.DSL().begin(); it != State.DSL().end(); ++it)
      setSizesBefore.push_back( (*it)->Size() );
    t_sync.Stop();
  }
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
  mprintf("    PARALLELANALYSIS: Will attempt to run current Analyses in parallel.\n\n"
          "*** THIS COMMAND IS STILL EXPERIMENTAL! ***\n\n");
  if (syncToMaster)
    mprintf("\tResulting data sets will be synced back to master.\n");
  // Naively divide up all analyses among threads.
  int my_start, my_stop;
  int nelts = Parallel::World().DivideAmongProcesses( my_start, my_stop, State.Analyses().size() );
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
  State.DFL().AllThreads_WriteAllDF();
  State.Analyses().Clear();
  if (syncToMaster) {
    t_sync.Start();
    // Check which sizes have changed.
    if (setSizesBefore.size() != State.DSL().size()) {
      mprintf("Warning: Number of sets have changed. Not attempting to sync sets to master.\n");
    } else {
      for (unsigned int idx = 0; idx != State.DSL().size(); idx++) {
        int setHasChanged = 0;
        if (!Parallel::World().Master()) {
          if (setSizesBefore[idx] != State.DSL()[idx]->Size()) {
            if (State.Debug() > 0)
              rprintf("Set '%s' size has changed from %u to %zu\n",
                      State.DSL()[idx]->legend(), setSizesBefore[idx], State.DSL()[idx]->Size());
            setHasChanged = 1;
          }
        }
        int totalChanged;
        Parallel::World().AllReduce(&totalChanged, &setHasChanged, 1, MPI_INT, MPI_SUM);
        if (totalChanged > 0) {
          if (totalChanged == 1) {
            int sourceRank = 0;
            if (setHasChanged == 1)
              setHasChanged = Parallel::World().Rank();
            Parallel::World().ReduceMaster(&sourceRank, &setHasChanged, 1, MPI_INT, MPI_SUM);
            if (State.Debug() > 0)
              mprintf("DEBUG: Need to sync '%s' from %i\n", State.DSL()[idx]->legend(), sourceRank);
            if (Parallel::World().Master())
              State.DSL()[idx]->RecvSet( sourceRank, Parallel::World() );
            else if (setHasChanged == Parallel::World().Rank())
              State.DSL()[idx]->SendSet( 0,          Parallel::World() );
          } else
            mprintf("Warning: '%s' exists on multiple threads. Not syncing.\n",
                    State.DSL()[idx]->legend());
        }
      }
    }
    t_sync.Stop();
  }
  t_total.Stop();
  if (syncToMaster)
    t_sync.WriteTiming(2, "Sync:", t_total.Total());
  t_total.WriteTiming(1, "Total:");
  return CpptrajState::OK;
}
#else /* MPI */
// =============================================================================
// NON-MPI CODE BELOW
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
#endif /* MPI */

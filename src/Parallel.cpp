#include <cstdio>
#include "Parallel.h"

/** MPI world communicator. */
Parallel::Comm Parallel::world_ = Parallel::Comm();

#ifdef MPI
// printMPIerr()
/** Wrapper for MPI_Error string.  */
void Parallel::printMPIerr(int err, const char *routineName, int rank) {
  int len, eclass;
  char buffer[1024];

  MPI_Error_string(err, buffer, &len);
  MPI_Error_class(err, &eclass);
  // Remove newlines from MPI error string
  for (int i = 0; i != len; i++)
    if (buffer[i] == '\n') buffer[i] = ':';
  fprintf(stderr, "[%i] MPI ERROR %d: %s: [%s]\n", rank, eclass, routineName, buffer);

  return;
}

// checkMPIerr()
/** \return 1 and print error message if MPI action failed. */
int Parallel::checkMPIerr(int err, const char *routineName, int rank) {
  if (err != MPI_SUCCESS) {
    printMPIerr(err, routineName, rank);
    return 1;
  }
  return 0;
}

// Parallel::Init()
int Parallel::Init(int argc, char** argv) {
  printf("Entering Parallel::Init\n");
  if ( MPI_Init(&argc, &argv) != MPI_SUCCESS ) {
    fprintf(stderr,"Error: Could not initialize MPI.\n");
    return 1;
  }
  world_ = Comm(MPI_COMM_WORLD);
//# ifdef PARALLEL_DEBUG_VERBOSE
//  parallel_debug_init();
//# endif
  return 0;
}

int Parallel::End() {
//# ifdef PARALLEL_DEBUG_VERBOSE
//  parallel_debug_end();
//# endif
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}

#else /* MPI */
// ----- NON-MPI VERSIONS OF ROUTINES ------------------------------------------
int Parallel::Init(int argc, char** argv) { return 0; }
int Parallel::End() { return 0; }
#endif

// ===== Parallel::Comm ========================================================
#ifdef MPI
/** CONSTRUCTOR: Use given MPI communicator */
Parallel::Comm::Comm(MPI_Comm commIn) : comm_(commIn), rank_(0), size_(0) {
  MPI_Comm_size(comm_, &size_);
  MPI_Comm_rank(comm_, &rank_);
}
#endif

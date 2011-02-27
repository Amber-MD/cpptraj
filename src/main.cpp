/*
 * CPPTRAJ: Rewrite of PTRAJ in C++
 * 2010 Daniel R. Roe
 */
#include "PtrajState.h"
#include "PtrajMpi.h"
#include "CpptrajStdio.h"
#ifndef CPPTRAJ_VERSION_STRING
#define CPPTRAJ_VERSION_STRING "V1.0.1"
#endif

void Usage(char *programName) {
  mprinterr("Usage: %s [-p Top1, -p Top2, ...] [-i Input] [-debug N]\n",programName);
  mprinterr("       %s Top1 Input\n",programName);
}

/*
 * ---=== CPPTRAJ MAIN ROUTINE ===---
 * Call parallel Init (does nothing if not a parallel run)
 * Process input from command line/inputfiles/stdin
 * Run
 */
int main(int argc, char **argv) {
  PtrajState State;
  int err;

  // Parallel Init: NOTE Should check for err
  parallel_init(argc,argv);

  mprintf("\nCPPTRAJ: Trajectory Analysis. %s\n",CPPTRAJ_VERSION_STRING);
  mprintf("    ___  ___  ___  ___\n");
  mprintf("     | \\/ | \\/ | \\/ | \n");
  mprintf("    _|_/\\_|_/\\_|_/\\_|_\n\n");
#ifdef MPI
  mprintf("Running on %i processors\n\n",worldsize);
#endif

  err = State.ProcessCmdLineArgs(argc,argv);
  switch ( err ) {
    case 0 : State.Run(); break;
    case 1 : Usage(argv[0]); break;
  }

  parallel_end();

  mprintf("\n");
  return 0;
}

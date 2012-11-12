/*
 * CPPTRAJ: Rewrite of PTRAJ in C++
 * 2010 Daniel R. Roe
 */
#include <cstdio>
#include "Cpptraj.h"
#include "MpiRoutines.h"
#ifndef CPPTRAJ_VERSION_STRING
#define CPPTRAJ_VERSION_STRING "V13.9.2b"
#define CPPTRAJ_INTERNAL_VERSION "V3.7.2b"
#endif

// ----------========== CPPTRAJ MAIN ROUTINE ==========----------
/// Main routine.
int main(int argc, char **argv) {
  Cpptraj State;
  printf("\nCPPTRAJ: Trajectory Analysis. %s\n",CPPTRAJ_VERSION_STRING);
  printf("    ___  ___  ___  ___\n");
  printf("     | \\/ | \\/ | \\/ | \n");
  printf("    _|_/\\_|_/\\_|_/\\_|_\n");
  // Parallel Init: NOTE Should check for err
  parallel_init(argc,argv);
#ifdef MPI
  if (worldrank==0) printf("Running on %i processors\n",worldsize);
#endif
  Cpptraj::Mode cmode = State.ProcessCmdLineArgs(argc,argv);
  switch ( cmode ) {
    case Cpptraj::C_OK          : State.Run(); break;
    case Cpptraj::C_INTERACTIVE : 
      if (State.Interactive() == Cpptraj::C_OK) State.Run(); 
      break;
    case Cpptraj::C_ERR         : 
    case Cpptraj::C_QUIT        : break;
  }

  parallel_end();

  printf("\n");
  return 0;
}


/*
 * CPPTRAJ: Rewrite of PTRAJ in C++
 * 2010 Daniel R. Roe
 */
//#include <cstdio>
#include "PtrajState.h"
#include "PtrajMpi.h"

void Usage(char *programName) {
  mprintf(stderr,"Usage: %s [-p Top1, -p Top2, ...] [-i Input] [-debug N]\n",programName);
  mprintf(stderr,"       %s Top1 Input\n",programName);
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

  mprintf(stdout,"\nCPPTRAJ: Trajectory Analysis.\n");
  mprintf(stdout,"    ___  ___  ___  ___\n");
  mprintf(stdout,"     | \\/ | \\/ | \\/ | \n");
  mprintf(stdout,"    _|_/\\_|_/\\_|_/\\_|_\n\n");
  // NOTE: Eventually allow stdin mode, like Ptraj
  //if (argc<2) {
  //  Usage(argv[0]);
  //  return 0;
  //}

  // Not set up for parallel yet. 
  //if (worldrank==0) {
    err = State.ProcessCmdLineArgs(argc,argv);
    switch ( err ) {
      case 0 : State.Run(); break;
      case 1 : Usage(argv[0]); break;
    }

  //}

  parallel_end();

  mprintf(stdout,"\n");
  return 0;
}

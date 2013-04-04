/*
 * CPPTRAJ: Rewrite of PTRAJ in C++
 * Copyright (c) 2010-2013 Daniel R. Roe
 * For license information see the LICENSE file.
 * For a full list of contributing authors see the README file.
 */
#include <unistd.h> // isatty
#include <cstdio>
#include "Cpptraj.h"
#include "MpiRoutines.h"
#include "Version.h"
// ----------========== CPPTRAJ MAIN ROUTINE ==========----------
/// Main routine.
int main(int argc, char **argv) {
  int err = 0;
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
    case Cpptraj::C_OK          : 
      err = State.Run(); break;
    case Cpptraj::C_INTERACTIVE :
      // Test if input is really from a console
      if ( isatty(fileno(stdin)) )
        cmode = State.Interactive();
      else
        cmode = State.ProcessInput(""); // "" means read from STDIN 
      if (cmode == Cpptraj::C_OK) 
        err = State.Run(); 
      else if (cmode == Cpptraj::C_ERR)
        err = 1;
      break;
    case Cpptraj::C_ERR         :
      err = 1;
    case Cpptraj::C_QUIT        : break;
  }
  parallel_end();
  printf("\n");
  return err;
}

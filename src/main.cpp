/*
 * CPPTRAJ: Rewrite of PTRAJ in C++
 * Copyright (c) 2010-2013 Daniel R. Roe
 * For license information see the LICENSE file.
 * For a full list of contributing authors see the README file.
 */
#include <unistd.h> // isatty
#include <cstdio> // stdin, fileno
#include "Cpptraj.h"
#include "MpiRoutines.h"
// ----------========== CPPTRAJ MAIN ROUTINE ==========----------
/// Main routine.
int main(int argc, char **argv) {
  int err = 0;
  Cpptraj State;
  // Parallel Init: NOTE Should check for err
  parallel_init(argc,argv);
  Cpptraj::Intro();
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
  if (worldrank==0) printf("\n");
  parallel_end();
  return err;
}

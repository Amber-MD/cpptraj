/*
 * CPPTRAJ: Rewrite of PTRAJ in C++
 * Copyright (c) 2010-2013 Daniel R. Roe
 * For license information see the LICENSE file.
 * For a full list of contributing authors see the README file.
 */
#include "Cpptraj.h"
#include "MpiRoutines.h"
// ----------========== CPPTRAJ MAIN ROUTINE ==========----------
/// Main routine.
int main(int argc, char **argv) {
  int err = 0;
  Cpptraj State;
  if (parallel_init(argc,argv) != 0) return 1;
  Cpptraj::Intro();
  Cpptraj::Mode cmode = State.ProcessCmdLineArgs(argc,argv);
  switch ( cmode ) {
    case Cpptraj::C_OK          :
      // If run has not yet been called, run now.
      if (State.Nrun() < 1) 
        err = State.Run(); 
      break;
    case Cpptraj::C_INTERACTIVE :
      if (State.Interactive() == Cpptraj::C_ERR)
        err = 1;
      break;
    case Cpptraj::C_ERR         :
      err = 1;
    case Cpptraj::C_QUIT        : break;
  }
  Cpptraj::Finalize(err);
  parallel_end();
  return err;
}

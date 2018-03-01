/*
 * CPPTRAJ: Rewrite of PTRAJ in C++
 * Copyright (c) 2010-2016 Daniel R. Roe
 * For license information see the LICENSE file.
 * For a full list of contributing authors see the README file.
 */
#ifdef FFTW_FFT
#include <fftw3.h>
#endif
#include "Cpptraj.h"
#include "Parallel.h"

/** The Cpptraj class lives here to ensure destructors are called
  * before Parallel::End();
  */
int CpptrajProgram(int argc, char** argv) {
  Cpptraj Program;
  int err = Program.RunCpptraj(argc, argv);
# ifdef FFTW_FFT
  // Ensure no debris from FFTW is left over.
  fftw_cleanup();
# endif
  return err;
}

// ----------========== CPPTRAJ MAIN ROUTINE ==========----------
/// Main routine.
int main(int argc, char **argv) {
  if (Parallel::Init(argc, argv) != 0) return 1;
  int err = CpptrajProgram(argc, argv);
  Parallel::End();
  return err;
}

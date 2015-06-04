/*
 * CPPTRAJ: Rewrite of PTRAJ in C++
 * Copyright (c) 2010-2014 Daniel R. Roe
 * For license information see the LICENSE file.
 * For a full list of contributing authors see the README file.
 */
#ifdef FFTW_FFT
#include <fftw3.h>
#endif
#include "Cpptraj.h"
#include "MpiRoutines.h"
// ----------========== CPPTRAJ MAIN ROUTINE ==========----------
/// Main routine.
int main(int argc, char **argv) {
  Cpptraj Program;
  if (parallel_init(argc,argv) != 0) return 1;
  int err = Program.RunCpptraj(argc, argv);
  parallel_end();
# ifdef FFTW_FFT
  // Ensure no debris from FFTW is left over.
  fftw_cleanup();
# endif
  return err;
}

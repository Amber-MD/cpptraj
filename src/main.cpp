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
#include "Parallel.h"
// ----------========== CPPTRAJ MAIN ROUTINE ==========----------
/// Main routine.
int main(int argc, char **argv) {
  if (Parallel::Init(argc,argv) != 0) return 1;
  Cpptraj Program;
  int err = Program.RunCpptraj(argc, argv);
# ifdef FFTW_FFT
  // Ensure no debris from FFTW is left over.
  fftw_cleanup();
# endif
  Parallel::End();
  return err;
}

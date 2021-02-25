// Test the RNGs in cpptraj with TestU01
#include <string.h>
#include <stdlib.h>
#include "TestU01.h"
#include "cpptraj_rng.h"

static void* rngptr = 0;

unsigned int test_rng(void)
{
  return num_cpptraj_rng( rngptr );
}

int argIs(const char* arg, const char* key)
{
  if (strcmp(arg, key) == 0) return 1;
  return 0;
}

int main(int argc, char** argv) {
  int iseed = 0;
  int icount = 0;
  int itype = 0;
  int imax = 10;
  // Parse options
  for (int iarg = 1; iarg < argc; iarg++) {
    if (argIs(argv[iarg], "-S")) {
      iseed = atoi( argv[++iarg] );
    } else if (argIs(argv[iarg], "-t")) {
      icount = atoi( argv[++iarg] );
    } else if (argIs(argv[iarg], "-r")) {
      itype = atoi( argv[++iarg] );
    } else if (argIs(argv[iarg], "--max")) {
      imax = atoi( argv[++iarg] );
    }/* else if (arg == "--mode") {
      std::string marg(argv[++iarg]);
      if (marg == "diehard")
        mode = DIEHARD;
      else if (marg == "shuffle")
        mode = SHUFFLE;
      else if (marg == "range")
        mode = RANGE;
      else {
        fprintf(stderr,"Error: Unrecognized mode: %s\n", marg.c_str());
        return 1;
      }
    }*/
  }

  // Setup RNG
  rngptr = get_cpptraj_rng( itype, iseed );
  
  //int it;
  //for (it = 0; it < icount; it++)
  //  printf("%10u\n", num_cpptraj_rng( rngptr ));

  // Create TestU01 PRNG object for our generator
  unif01_Gen* gen = unif01_CreateExternGenBits((char*)"CPPTRAJ", test_rng);

  // Run the tests.
  bbattery_SmallCrush(gen);

  // Clean up.
  unif01_DeleteExternGenBits(gen);

  destroy_cpptraj_rng( rngptr );
  return 0;
} 

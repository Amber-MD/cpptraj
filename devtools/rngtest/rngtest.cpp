#include <limits>
#include <string>
#include <cstdio>
#include <cstdlib>
#include "../../src/Random.h"
#include "../../src/CpptrajStdio.h"

// Test the random number generators in cpptraj
int main(int argc, char** argv) {
  SetWorldSilent(true);
  Random_Number RNG;

  int iseed = 0;
  int icount = 0;
  int itype = 0;
  // Parse options
  for (int iarg = 1; iarg < argc; iarg++) {
    std::string arg( argv[iarg] );
    if (arg == "-S") {
      iseed = atoi( argv[++iarg] );
    } else if (arg == "-t") {
      icount = atoi( argv[++iarg] );
    } else if (arg == "-r") {
      itype = atoi( argv[++iarg] );
    }
  }
  printf("#==================================================================\n");
  printf("# Cpptraj generator %i  seed = %i                                  \n", itype, iseed);
  printf("#==================================================================\n");
  printf("type: d\n");
  printf("count: %i\n", icount);
  printf("numbit: 32\n");
  //printf("\tSeed  : %i\n", iseed);
  //printf("\tCount : %i\n", icount);
  //printf("\tType  : %i\n", itype);

  Random_Number::RngType rt = (Random_Number::RngType)itype;
  Random_Number::SetDefaultRng(rt);
  RNG.rn_set( iseed );

  double d_imax = (double)std::numeric_limits<int>::max();
  for (int ic = 0; ic < icount; ic++) {
    double rn = RNG.rn_gen() * d_imax;
    int in = (int)rn;
    printf("%10i\n", in);
  }

  printf("\n");

  return 0;
}

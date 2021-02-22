#include <limits>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include "../../src/Random.h"
#include "../../src/CpptrajStdio.h"

static int DieHard_Tests(Random_Number const& RNG, int itype, int iseed, int icount) {
  printf("#==================================================================\n");
  printf("# Cpptraj generator %s  seed = %i                                  \n",
         Random_Number::CurrentDefaultRngStr(), iseed);
  printf("#==================================================================\n");
  printf("type: d\n");
  printf("count: %i\n", icount);
  printf("numbit: 32\n");

/*
  double d_imax = (double)std::numeric_limits<int>::max();
  for (int ic = 0; ic < icount; ic++) {
    double rn = RNG.rn_gen() * d_imax;
    int in = (int)rn;
    printf("%10i\n", in);
  }*/
  for (int ic = 0; ic < icount; ic++) {
    printf("%10u\n", RNG.rn_num());
  }

  //printf("\n");
  return 0;
}

static int Shuffle_Tests(Random_Number const& RNG, int itype, int iseed, int icount) {
  printf("Shuffle tests\n");
  printf("\tSeed  : %i\n", iseed);
  printf("\tCount : %i\n", icount);
  printf("\tType  : %i\n", itype);

  std::vector<int> nums;
  nums.reserve(icount);
  for (int i = 0; i < icount; i++)
    nums.push_back( i );

  RNG.ShufflePoints( nums );
  printf("Shuffled points:");
  for (std::vector<int>::const_iterator it = nums.begin(); it != nums.end(); ++it)
    printf(" %i", *it);
  printf("\n");
  return 0;
}

// Test the random number generators in cpptraj
int main(int argc, char** argv) {
  SetWorldSilent(true);
  Random_Number RNG;

  enum ModeType {DIEHARD = 0, SHUFFLE };
  ModeType mode = DIEHARD;
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
    } else if (arg == "--mode") {
      std::string marg(argv[++iarg]);
      if (marg == "diehard")
        mode = DIEHARD;
      else if (marg == "shuffle")
        mode = SHUFFLE;
      else {
        fprintf(stderr,"Error: Unrecognized mode: %s\n", marg.c_str());
        return 1;
      }
    }
  }

  // Setup RNG
  Random_Number::RngType rt = (Random_Number::RngType)itype;
  Random_Number::SetDefaultRng(rt);
  if (RNG.rn_set( iseed )) {
    fprintf(stderr, "Error: Could not set up RNG.\n");
    return 1;
  }
  int err = 0;

  if (mode == DIEHARD) {
    err = DieHard_Tests(RNG, itype, iseed, icount);
  } else if (mode == SHUFFLE) {
    err = Shuffle_Tests(RNG, itype, iseed, icount);
  }

  return err;
}

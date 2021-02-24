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
  printf("Shuffle tests (%s)\n", Random_Number::CurrentDefaultRngStr());
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

static int Range_Tests(Random_Number const& RNG, int itype, int iseed, int icount) {
  printf("Range tests (%s)\n", Random_Number::CurrentDefaultRngStr());
  printf("\tSeed  : %i\n", iseed);
  printf("\tCount : %i\n", icount);
  printf("\tType  : %i\n", itype);

  unsigned int min = 1;
  unsigned int max = 10;
  std::vector<int> Counts( max-min+1, 0 );
  unsigned int oob = 0;
  for (int i = 0; i < icount; i++) {
    unsigned int rn = RNG.rn_num_interval(min, max);
    long int idx = (long int)rn - (long int)min;
    if (idx < 0 || idx >= (long int)Counts.size()) {
      printf("Random # %u out of bounds (%li)\n", rn, idx);
      oob++;
    } else {
      Counts[idx]++;
    }
  }

  printf("Final counts (%u out of bounds):\n", oob);
  unsigned int num = min;
  for (std::vector<int>::const_iterator it = Counts.begin(); it != Counts.end(); ++it, ++num)
    printf("\t%10u : %10u\n", num, *it);
  return 0;
}

// Test the random number generators in cpptraj
int main(int argc, char** argv) {
  SetWorldSilent(true);
  Random_Number RNG;

  enum ModeType {DIEHARD = 0, SHUFFLE, RANGE };
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
      else if (marg == "range")
        mode = RANGE;
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
  } else if (mode == RANGE) {
    err = Range_Tests(RNG, itype, iseed, icount);
  }

  return err;
}

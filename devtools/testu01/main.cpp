// Test the RNGs in cpptraj with TestU01
#include <string>
#include <cstdlib>
#include "TestU01.h"
#include "../../src/Random.h"

static Random_Number* rng_ = 0;

unsigned int test_rng(void)
{
  return rng_->rn_num();
}

int main(int argc, char** argv) {
  int iseed = 0;
  int icount = 0;
  int itype = 0;
  int imax = 10;
  // Parse options
  for (int iarg = 1; iarg < argc; iarg++) {
    std::string arg( argv[iarg] );
    if (arg == "-S") {
      iseed = atoi( argv[++iarg] );
    } else if (arg == "-t") {
      icount = atoi( argv[++iarg] );
    } else if (arg == "-r") {
      itype = atoi( argv[++iarg] );
    } else if (arg == "--max") {
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
  Random_Number::RngType rt = (Random_Number::RngType)itype;
  Random_Number::SetDefaultRng(rt);
  rng_ = new Random_Number();
  if (rng_->rn_set( iseed )) {
    fprintf(stderr, "Error: Could not set up RNG.\n");
    return 1;
  }

  if (rng_ != 0) delete rng_;
  return 0;
} 

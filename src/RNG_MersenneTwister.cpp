#include "RNG_MersenneTwister.h"

double Cpptraj::RNG_MersenneTwister::Generate() {
# ifdef C11_SUPPORT
  unsigned urn = gen_();
  return (double)urn / (double)gen_.max();
# else
  return 0;
# endif
}

int Cpptraj::RNG_MersenneTwister::SetupRng() {
# ifdef C11_SUPPORT
  gen_.seed( Seed() );
  return 0;
# else
  mprinterr("Error: Cpptraj was compiled without C++11 support. Cannot use Mersenne twister.\n");
  return 1;
# endif
}

#include "RNG_PCG32.h"
#include "CpptrajStdio.h"

int Cpptraj::RNG_PCG32::SetupRng() {
# ifdef _MSC_VER
  mprinterr("Error: PCG32 generator does not yet work under Visual Studio\n");
  return 1;
# else
  rng_ = pcg32( Seed() );
  return 0;
# endif
}
/*
double Cpptraj::RNG_PCG32::Generate() {
# ifdef _MSC_VER
  return 0;
# else
  unsigned int rn = rng_();
  return (double)rn / (double)rng_.max();
# endif
}*/

unsigned int Cpptraj::RNG_PCG32::Number() {
# ifdef _MSC_VER
  return 0;
# else
  return rng_();
# endif
}

#include "RNG_Stdlib.h"
#include <cstdlib>

Cpptraj::RNG_Stdlib::RNG_Stdlib() {}

int Cpptraj::RNG_Stdlib::SetupRng() {
  srand( Seed() );
  return 0;
}

double Cpptraj::RNG_Stdlib::Generate() {
  const double drand_max = (double)RAND_MAX;
  int rn = rand();
  double drn = (double)rn;

  return (drn / drand_max);
}

unsigned int Cpptraj::RNG_Stdlib::Number() {
  return rand();
}

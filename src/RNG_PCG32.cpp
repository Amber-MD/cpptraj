#include "RNG_PCG32.h"

int Cpptraj::RNG_PCG32::SetupRng() {
  rng_ = pcg32( Seed() );
}

double Cpptraj::RNG_PCG32::Generate() {
  unsigned int rn = rng_();
  return (double)rn / (double)rng_.max();
}

unsigned int Cpptraj::RNG_PCG32::Number() {
  return rng_();
}

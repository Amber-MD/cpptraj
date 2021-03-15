#include "xoshiro128plusplus.h"
#include "RNG_Xoshiro128pp.h"
#include <limits>

int Cpptraj::RNG_Xoshiro128pp::SetupRng() {
  x128pp_seed( Seed() );
  return 0;
}

unsigned int Cpptraj::RNG_Xoshiro128pp::Number() {
  return x128pp_next();
}
/*
double Cpptraj::RNG_Xoshiro128pp::Generate() {
  unsigned int rn = Number();
  return (double)rn / (double)std::numeric_limits<unsigned int>::max();
}*/

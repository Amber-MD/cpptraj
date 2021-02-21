#include "CpptrajStdio.h"
#include "Random.h"
#include "RNG.h"
#include "RNG_Marsaglia.h"
#include "RNG_Stdlib.h"
#include "RNG_MersenneTwister.h"
#include "RNG_PCG32.h"
#include "RNG_Xoshiro128pp.h"

Random_Number::RngType Random_Number::defaultType_ = MARSAGLIAS;

int Random_Number::defaultSeed_ = 71277; // AMBER default seed

void Random_Number::SetDefaultRng(RngType r) {
  defaultType_ = r;
}

void Random_Number::SetDefaultSeed(int i) {
  defaultSeed_ = i;
}

/** CONSTRUCTOR */
Random_Number::Random_Number() :
  rng_(0)
{}

/** DESTRUCTOR */
Random_Number::~Random_Number() {
  if (rng_ != 0) delete rng_;
}

/** Allocate RNG. */
void Random_Number::allocateRng() {
  if (rng_ != 0) delete rng_;
  rng_ = 0;
  switch (defaultType_) {
    case MARSAGLIAS :
      mprintf("\tRNG: Marsaglia\n");
      rng_ = new Cpptraj::RNG_Marsaglia();
      break;
    case STDLIB     :
      mprintf("\tRNG: C stdlib\n");
      rng_ = new Cpptraj::RNG_Stdlib();
      break;
    case MERSENNE_TWISTER :
      mprintf("\tRNG: Mersenne twister\n");
      rng_ = new Cpptraj::RNG_MersenneTwister();
      break;
    case PCG32            :
      mprintf("\tRNG: PCG 32\n");
      rng_ = new Cpptraj::RNG_PCG32();
      break;
    case XOSHIRO128PP     :
      mprintf("\tRNG: Xoshiro128++\n");
      rng_ = new Cpptraj::RNG_Xoshiro128pp();
    break;
  } 
}

/** Initialize RNG. */
int Random_Number::rn_set(int seedIn) {
  allocateRng();
  return rng_->Set_Seed( seedIn );
}

/** Initialize with default seed. */
//void Random_Number::rn_set() {
//  rng_->Set_Seed();
//}

/** Generate random number. */
double Random_Number::rn_gen() const {
  return rng_->Generate();
}

/** Generate a random integer. */
int Random_Number::rn_num() const {
  return rng_->Number();
}

/** Generate pseudo-random Gaussian sequence. */
double Random_Number::rn_gauss(double am, double sd) const {
  return rng_->GenerateGauss(am, sd);
}

/** \return true if RNG has been set up. */
bool Random_Number::IsSet() const {
  return rng_->IsSet();
}

/** \return Value of RNG seed. */
int Random_Number::Seed() const {
  return rng_->Seed();
}

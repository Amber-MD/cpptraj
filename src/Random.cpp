#include <cmath> //sqrt, log
#include "CpptrajStdio.h"
#include "Random.h"
#include "RNG.h"
#include "RNG_Marsaglia.h"
#include "RNG_Stdlib.h"
#include "RNG_MersenneTwister.h"
#include "RNG_PCG32.h"
#include "RNG_Xoshiro128pp.h"

/** Starting default type. */
Random_Number::RngType Random_Number::defaultType_ = Random_Number::MARSAGLIA;
// TODO Test this generator, then enable. Should be a better RNG than Marsaglia
//Random_Number::RngType Random_Number::defaultType_ = Random_Number::XOSHIRO128PP;

/** Starting default seed. */
//int Random_Number::defaultSeed_ = 71277; // AMBER default seed

/** Set the default RNG type. */
void Random_Number::SetDefaultRng(RngType r) {
  defaultType_ = r;
}

/** Set the default seed. */
//void Random_Number::SetDefaultSeed(int i) {
//  defaultSeed_ = i;
//}

/** CONSTRUCTOR */
Random_Number::Random_Number() :
  rng_(0)
{}

/** DESTRUCTOR */
Random_Number::~Random_Number() {
  if (rng_ != 0) delete rng_;
}

/** \return Description of current default RNG. */
const char* Random_Number::CurrentDefaultRngStr() {
  const char* str = 0;
  switch (defaultType_) {
    case MARSAGLIA        : str = "Marsaglia"; break;
    case STDLIB           : str = "C stdlib"; break;
    case MERSENNE_TWISTER : str = "Mersenne Twister (mt19937)"; break;
    case PCG32            : str = "Permuted Congruential Generator (32 bit)"; break;
    case XOSHIRO128PP     : str = "Xoshiro128++"; break;
  }
  return str;
}

/** Allocate RNG according to the current default type. */
void Random_Number::allocateRng() {
  if (rng_ != 0) delete rng_;
  rng_ = 0;
  switch (defaultType_) {
    case MARSAGLIA  :
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

/** Allocate RNG, initialize with given seed. */
int Random_Number::rn_set(int seedIn) {
  allocateRng();
  return rng_->Set_Seed( seedIn );
}

/** Initialize with default seed. */
//void Random_Number::rn_set() {
//  rng_->Set_Seed();
//}

/** Generate random number between 0 and 1. */
double Random_Number::rn_gen() const {
  return rng_->Generate();
}

/** Generate a random integer. */
unsigned int Random_Number::rn_num() const {
  return rng_->Number();
}

/** Generate a random integer on interval imin to imax. 
  * NOTE: Unlike rn_num(), this returns an integer to allow intervals of
  *       negative numbers.
  */
int Random_Number::rn_num_interval_signed(int imin, int imax) const {
  int iwidth = imax - imin;
  if (iwidth < 1) {
    mprinterr("Internal Error: Random_Number::rn_num_interval_signed(): Interval width <= 0.\n");
    return 0;
  }
  unsigned int uwidth = (unsigned int)iwidth + 1;
  //unsigned int umod = (rn_num() % uwidth);
  unsigned int umod = rng_->Number_UpTo( uwidth );
  return (int)umod + imin;
}

/** Generate a random integer on interval umin to umax. */
unsigned int Random_Number::rn_num_interval(unsigned int umin, unsigned int umax) const {
  if (umax > umin) {
   unsigned int uwidth = (umax - umin) + 1;
   //unsigned int umod = (rn_num() % uwidth);
   unsigned int umod = rng_->Number_UpTo( uwidth );
   return umod + umin; 
  } else {
    mprinterr("Internal Error: Random_Number::rn_num_interval(): Interval max <= interval min.\n");
  }
  return 0;
}

/** This uses Generate() to generate random numbers between 0 and 1
  * in a Gaussian distribution, with mean "am" and standard deviation "sd".
  */
double Random_Number::GenerateGauss(double am, double sd) const {
  if (!IsSet()) {
    mprinterr("Error: random number generator not initialized.");
    return -1.0;
  }
  // Use the method of Box and Muller.
  // For some applications, one could use both "v" and "veven" in random
  // sequence; but this won't work for most things we need (e.g. Langevin
  // dynamics,) since the two adjacent variables are themselves highly
  // correlated.  Hence we will just use the first ("v") variable.

  // Get two random numbers, even on (-1,1):
  bool generate = true;
  while (generate) {
    double uni = rng_->Generate();
    double zeta1 = uni + uni - 1.0;

    uni = rng_->Generate();
    double zeta2 = uni + uni - 1.0;

    double tmp1 = zeta1 * zeta1 + zeta2 * zeta2;

    if (tmp1 < 1.0 && tmp1 != 0.0) {
        double tmp2 = sd * sqrt(-2.0 * log(tmp1)/tmp1);
        return (zeta1 * tmp2 + am);
    }
  }
  return 0.0;
}

/** Generate pseudo-random Gaussian sequence. */ // TODO deprecate this version
double Random_Number::rn_gauss(double am, double sd) const {
  return GenerateGauss(am, sd);
}

/** Use modern version of the Fisher-Yates shuffle to randomly reorder the
  * given points.
  */
void Random_Number::ShufflePoints( std::vector<int>& PointIndices ) const {
  for (unsigned int i = PointIndices.size() - 1; i != 1; i--)
  { // 0 <= j <= i
    unsigned int j = rn_num_interval(0, i);
    int temp = PointIndices[j];
    PointIndices[j] = PointIndices[i];
    PointIndices[i] = temp;
  }
}

/** \return true if RNG has been set up. */
bool Random_Number::IsSet() const {
  return rng_->IsSet();
}

/** \return Value of RNG seed. */
int Random_Number::Seed() const {
  return rng_->Seed();
}

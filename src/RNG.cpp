#include <ctime> // clock
#include "RNG.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
Cpptraj::RNG::RNG() : iseed_(-1) {}

/** Set up the seed, then perform any specific setup required by the RNG. */
int Cpptraj::RNG::Set_Seed(int iseedIn) {
  // If iseed <= 0 set the seed based on wallclock time
  if (iseedIn <= 0 ) {
    clock_t wallclock = clock();
    iseed_ = (int)wallclock;
    mprintf("Random_Number: seed is <= 0, using wallclock time as seed (%i)\n", iseed_);
  } else
    iseed_ = iseedIn;
  // Set up specific to RNG
  if (SetupRng()) {
    mprinterr("Error: RNG setup failed.\n");
    return 1;
  }
  return 0;
}

/** Generate random number from 0 to a given max. Use multiply and 
 *  bit-shift to try to avoid modulo bias.
 */
unsigned int Cpptraj::RNG::Number_UpTo(unsigned int exclusiveMax)
{
  if (exclusiveMax == 0) {
    mprinterr("Internal Error: RandomNum::rn_num_max(): max is zero.\n");
    return 0;
  }
  return (unsigned int)(((unsigned long)exclusiveMax * Number()) >> 32);
}

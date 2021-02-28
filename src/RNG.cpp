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

/** Generate random number from 0 to a given max. Try to avoid modulo bias.
 */
unsigned int Cpptraj::RNG::Number_UpTo(unsigned int exclusiveMax)
{
  if (exclusiveMax <= 1u<<12) {
    // Range is up to 2^12
    if (exclusiveMax == 0) {
      mprinterr("Internal Error: RandomNum::rn_num_max(): max is zero.\n");
      return 0;
    }
    // Use multiply and bit-shift.
    return (unsigned int)(((unsigned long)exclusiveMax * Number()) >> 32);
  } else {
    // Make limit = 2^32 - N (i.e. exclusiveMax), compare r-(r%N) to limit.
    // If less than the limit it is a multiple of N and should be free of
    // modulo bias.
    unsigned int r, v, limit = (unsigned int)-(int)exclusiveMax;
    do
    {
      r = Number();
      v = r % exclusiveMax;
    } while(r-v > limit);
    return v;
  }
}

#define M_RAN_INVM32 2.32830643653869628906e-010
/** Generate a random number between 0 and 1 from 32bit unsigned int.
  * Uses code from Doornik, ACM Transactions on Mathematical Software, 2006.
  */
double Cpptraj::RNG::Generate() {
  unsigned int uiRan = Number();
  return (int)uiRan * M_RAN_INVM32 + 0.5;
}

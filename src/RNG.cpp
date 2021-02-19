#include <ctime> // clock
#include <cmath> // sqrt, log
#include "RNG.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
Cpptraj::RNG::RNG() : iseed_(-1) {}

/** Set up the seed, then perform any specific setup required by the RNG. */
void Cpptraj::RNG::rn_set(int iseedIn) {
  // If iseed <= 0 set the seed based on wallclock time
  if (iseedIn <= 0 ) {
    clock_t wallclock = clock();
    iseed_ = (int)wallclock;
    mprintf("Random_Number: seed is <= 0, using wallclock time as seed (%i)\n", iseed_);
  } else
    iseed_ = iseedIn;

  if (SetupRng()) {
    mprinterr("Error: RNG setup failed.\n");
  }
}

/** This is a version of Generate() that adds the constraint of a Gaussian
  * distribution, with mean "am" and standard deviation "sd". This 
  * requires rn_set() to have been called first, and "uses up" the
  * same sequence that Generate() does.
  */
double Cpptraj::RNG::rn_gauss(double am, double sd) {
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
    double uni = Generate();
    double zeta1 = uni + uni - 1.0;

    uni = Generate();
    double zeta2 = uni + uni - 1.0;

    double tmp1 = zeta1 * zeta1 + zeta2 * zeta2;

    if (tmp1 < 1.0 && tmp1 != 0.0) {
        double tmp2 = sd * sqrt(-2.0 * log(tmp1)/tmp1);
        return (zeta1 * tmp2 + am);
    }
  }
  return 0.0;
}

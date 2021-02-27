#include <cstdlib> // div
#include <limits>
#include "RNG_Marsaglia.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
/// Only sets setup state to false. Actual setup is rn_set
Cpptraj::RNG_Marsaglia::RNG_Marsaglia() {
  RN_generator.iseed = -1;
}

/** Max value for this RNG, 2^24 bits */
const double Cpptraj::RNG_Marsaglia::rng_max_ = 16777216.0;

/** Initialization routine for Marsaglias random number generator
  * as implemented in Amber 3.0 Rev A by George Seibel. See doc in amrand.
  *
  * Testing:  Call amrset with iseed = 54185253.  This should result
  * in is1 = 1802 and is2 = 9373.  Call amrand 20000 times, then six
  * more times, printing the six random numbers * 2**24 (ie, 4096*4096)
  * They should be: (6f12.1)
  * 6533892.0  14220222.0  7275067.0  6172232.0  8354498.0  10633180.0
  */
int Cpptraj::RNG_Marsaglia::SetupRng() {
  RN_generator.iseed = Seed();
  // Two internal seeds used in Marsaglia algorithm:
  int is1, is2;
  const int is1max = 31328;
  const int is2max = 30081;

  div_t divresult;  

  // Construct two internal seeds from single unbound seed.
  // is1 and is2 are quotient and remainder of iseed/IS2MAX.
  // Add one to keep zero and one results from both mapping to one.
  divresult = div( RN_generator.iseed, is2max );
  divresult.quot++;
  divresult.rem++;
  // Ensure range: 1 <= is1 <= is1max
  if ( divresult.quot > 1 )
    is1 = divresult.quot;
  else
    is1 = 1;
  if ( is1max < is1 )
    is1 = is1max;
  // Ensure range: 1 <= is2 <= is2max
  if ( divresult.rem > 1 )
    is2 = divresult.rem;
  else
    is2 = 1;
  if ( is2max < is2 )
    is2 = is2max;

  int i = ((is1/177) % 177) + 2;
  int j = (is1       % 177) + 2;
  int k = ((is2/169) % 178) + 1;
  int l = (is2       % 169);
    
  for (int ii = 0; ii < 97; ii++) { 
    double s = 0.0;
    double t = 0.5;
    for (int jj = 0; jj < 24; jj++) {
      //m = mod(mod(i*j, 179)*k, 179)
      int mtemp = (i*j) % 179;
      int m = (mtemp*k) % 179; 
      i = j;
      j = k;
      k = m;
      l = (((53*l) + 1) % 169);
      //if (mod(l*m, 64) .ge. 32) s = s + t
      if ( ((l*m) % 64) >= 32 ) s += t;
      //t = 0.5d0 * t
      t *= 0.5;
    }
    RN_generator.u[ii] = s;
  }

  RN_generator.c  = 362436.0   / rng_max_;
  RN_generator.cd = 7654321.0  / rng_max_;
  RN_generator.cm = 16777213.0 / rng_max_;

  RN_generator.i97 = 96;
  RN_generator.j97 = 32;
  return 0;
}

/** \author Portable Random number generator by George Marsaglia
  * Original author: Amber 3.0 Rev A implementation by George Seibel
  *
  * This random number generator originally appeared in 'Toward a Universal
  * Random Number Generator' by George Marsaglia and Arif Zaman.  Florida
  * State University Report: FSU-SCRI-87-50 (1987)
  *
  * It was later modified by F. James and published in 'A Review of Pseudo-
  * random Number Generators'
  *
  * This is claimed to be the best known random number generator available.
  * It passes ALL of the tests for random number generators and has a
  * period of 2^144, is completely portable (gives bit identical results on
  * all machines with at least 24-bit mantissas in the floating point
  * representation).
  *
  * The algorithm is a combination of a Fibonacci sequence (with lags of 97
  * and 33, and operation "subtraction plus one, modulo one") and an
  * "arithmetic sequence" (using subtraction).
  *
  * \return A random number between 0.0 and 1.0
  * \return -1.0 if the random number generator is not initialized.
  */
double Cpptraj::RNG_Marsaglia::Generate() {
  if (!IsSet()) { 
    mprinterr("Error: random number generator not initialized.");
    return -1.0;
  }
  
  double uni = RN_generator.u[RN_generator.i97] - RN_generator.u[RN_generator.j97];
  if (uni < 0.0) uni += 1.0;
  RN_generator.u[RN_generator.i97] = uni;
  RN_generator.i97--;
  if (RN_generator.i97 == -1) RN_generator.i97 = 96;
  RN_generator.j97--;
  if (RN_generator.j97 == -1) RN_generator.j97 = 96;
  RN_generator.c -= RN_generator.cd;
  if (RN_generator.c < 0.0) RN_generator.c += RN_generator.cm;
  uni -= RN_generator.c;
  if (uni < 0.0) uni += 1.0; 
  return uni;
}

unsigned int Cpptraj::RNG_Marsaglia::Number() {
  double rn = Generate();
  return (unsigned int)(rn * rng_max_);
}

unsigned int Cpptraj::RNG_Marsaglia::Number_UpTo(unsigned int exclusiveMax) {
  return (unsigned int)(Generate() * (double)exclusiveMax);
}

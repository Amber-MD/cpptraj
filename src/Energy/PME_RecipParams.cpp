#include "PME_RecipParams.h"
#include "../Box.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"

using namespace Cpptraj::Energy;

/** CONSTRUCTOR - NFFT=-1 means auto-assign. */
PME_RecipParams::PME_RecipParams() :
  order_(6),
  debug_(0)
{
  nfft_[0] = -1;
  nfft_[1] = -1;
  nfft_[2] = -1;
}

/** Initialize PME recip options. */
int PME_RecipParams::InitRecip(EwaldOptions const& pmeOpts, int debugIn) {
  debug_ = debugIn;
  nfft_[0] = pmeOpts.Nfft1();
  nfft_[1] = pmeOpts.Nfft2();
  nfft_[2] = pmeOpts.Nfft3();
  order_ = pmeOpts.SplineOrder();
  // Set defaults if necessary
  if (order_ < 1) order_ = 6;
  return 0;
}

/** Print options to stdout. */
void PME_RecipParams::PrintRecipOpts() const {
  //mprintf("\tRecip opts (distance kernel exponent= %i\n", distKernelExponent_);
  mprintf("\t  Bspline order= %i\n", order_);
  mprintf("\t ");
  for (int i = 0; i != 3; i++)
    if (nfft_[i] == -1)
      mprintf(" NFFT%i=auto", i+1);
    else
      mprintf(" NFFT%i=%i", i+1, nfft_[i]);
  mprintf("\n");
}

/** \return true if given number is a product of powers of 2, 3, or 5. */
bool PME_RecipParams::check_prime_factors(int nIn) {
  if (nIn == 1) return true;
  int NL = nIn;
  int NQ;
  // First divide down by 2
  while (NL > 0) {
    NQ = NL / 2;
    if (NQ * 2 != NL) break;
    if (NQ == 1) return true;
    NL = NQ;
  }
  // Next try 3
  while (NL > 0) {
    NQ = NL / 3;
    if (NQ * 3 != NL) break;
    if (NQ == 1) return true;
    NL = NQ;
  }
  // Last try 5
  while (NL > 0) {
    NQ = NL / 5;
    if (NQ * 5 != NL) break;
    if (NQ == 1) return true;
    NL = NQ;
  }
  return false;
}

/** Compute the ceiling of len that is also a product of powers of 2, 3, 5.
  * Use check_prime_factors to get the smallest integer greater or equal
  * than len which is decomposable into powers of 2, 3, 5.
  */
int PME_RecipParams::ComputeNFFT(double len) {
  int mval = (int)len - 1;
  for (int i = 0; i < 100; i++) {
    mval += 1;
    // Sanity check
    if (mval < 1) {
      mprinterr("Error: Bad box length %g, cannot get NFFT value.\n", len);
      return 0;
    }
    if (check_prime_factors(mval))
      return mval;
  }
  mprinterr("Error: Failed to get good FFT array size for length %g Ang.\n", len);
  return 0;
}

/** Given a box, determine number of FFT grid points in each dimension. */
int PME_RecipParams::DetermineNfft(int& nfft1, int& nfft2, int& nfft3, Box const& boxIn) const
{
  // FIXME Should this be checking if NFFT would change instead?
  nfft1 = nfft_[0];
  nfft2 = nfft_[1];
  nfft3 = nfft_[2];
  if (nfft1 < 1) {
    // Need even dimension for X direction
    nfft1 = ComputeNFFT( (boxIn.Param(Box::X) + 1.0) * 0.5 );
    nfft1 *= 2;
  }
  if (nfft2 < 1)
    nfft2 = ComputeNFFT( boxIn.Param(Box::Y) );
  if (nfft3 < 1)
    nfft3 = ComputeNFFT( boxIn.Param(Box::Z) );

  if (nfft1 < 1 || nfft2 < 1 || nfft3 < 1) {
    mprinterr("Error: Bad NFFT values: %i %i %i\n", nfft1, nfft2, nfft3);
    return 1;
  }
  if (debug_ > 0) mprintf("DEBUG: NFFTs: %i %i %i\n", nfft1, nfft2, nfft3);

  return 0;
}



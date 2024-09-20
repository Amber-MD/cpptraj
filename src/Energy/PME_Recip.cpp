#include "PME_Recip.h"
#include "../AtomMask.h"
#include "../Box.h"
#include "../CpptrajStdio.h"
#include "../Frame.h"

typedef helpme::Matrix<double> Mat;

using namespace Cpptraj::Energy;

PME_Recip::PME_Recip(Type typeIn) :
  debug_(0)
{
  switch (typeIn) {
    case COULOMB:
      distKernelExponent_ = 1;
      scaleFac_ = 1.0;
      break;
    case LJ:
      distKernelExponent_ = 6;
      scaleFac_ = -1.0;
      break;
  }
}

void PME_Recip::SetDebug(int d) { debug_ = d; }

/** \return true if given number is a product of powers of 2, 3, or 5. */
bool PME_Recip::check_prime_factors(int nIn) {
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
int PME_Recip::ComputeNFFT(double len) {
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
int PME_Recip::DetermineNfft(int& nfft1, int& nfft2, int& nfft3, Box const& boxIn) const
{
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

/** \return Reciprocal space part of PME energy calc. */
// TODO currently helPME needs the coords/charge arrays to be non-const, need to fix that
double PME_Recip::Recip_ParticleMesh(Darray& coordsDin, Box const& boxIn, Darray& ChargeIn,
                                     const int* nfftIn, double ew_coeffIn, int orderIn)
{
  t_recip_.Start();
  // This essentially makes coordsD and chargesD point to arrays.
  Mat coordsD(&coordsDin[0], ChargeIn.size(), 3);
  Mat chargesD(&ChargeIn[0], ChargeIn.size(), 1);
  int nfft1 = -1;
  int nfft2 = -1;
  int nfft3 = -1;
  if (nfftIn != 0) {
    nfft1 = nfftIn[0];
    nfft2 = nfftIn[1];
    nfft3 = nfftIn[2];
  }
  if ( DetermineNfft(nfft1, nfft2, nfft3, boxIn) ) {
    mprinterr("Error: Could not determine FFT grid spacing.\n");
    return 0.0;
  }
  // Instantiate double precision PME object
  // Args: 1 = Exponent of the distance kernel: 1 for Coulomb
  //       2 = Kappa
  //       3 = Spline order
  //       4 = nfft1
  //       5 = nfft2
  //       6 = nfft3
  //       7 = scale factor to be applied to all computed energies and derivatives thereof
  //       8 = max # threads to use for each MPI instance; 0 = all available threads used.
  // NOTE: Scale factor for Charmm is 332.0716
  // NOTE: The electrostatic constant has been baked into the Charge_ array already.
  //auto pme_object = std::unique_ptr<PMEInstanceD>(new PMEInstanceD());
  pme_object_.setup(distKernelExponent_, ew_coeffIn, orderIn, nfft1, nfft2, nfft3, scaleFac_, 0);
  // Sets the unit cell lattice vectors, with units consistent with those used to specify coordinates.
  // Args: 1 = the A lattice parameter in units consistent with the coordinates.
  //       2 = the B lattice parameter in units consistent with the coordinates.
  //       3 = the C lattice parameter in units consistent with the coordinates.
  //       4 = the alpha lattice parameter in degrees.
  //       5 = the beta lattice parameter in degrees.
  //       6 = the gamma lattice parameter in degrees.
  //       7 = lattice type
  pme_object_.setLatticeVectors(boxIn.Param(Box::X), boxIn.Param(Box::Y), boxIn.Param(Box::Z),
                                boxIn.Param(Box::ALPHA), boxIn.Param(Box::BETA), boxIn.Param(Box::GAMMA),
                                PMEInstanceD::LatticeType::XAligned);
  double erecip = pme_object_.computeERec(0, chargesD, coordsD);

  t_recip_.Stop();
  return erecip;
}


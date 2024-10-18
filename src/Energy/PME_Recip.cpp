#include "PME_Recip.h"
#ifdef LIBPME
#include "../Box.h"
#include "../CpptrajStdio.h"

typedef helpme::Matrix<double> Mat;

using namespace Cpptraj::Energy;

PME_Recip::PME_Recip(Type typeIn)
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

/** Init */
int PME_Recip::InitRecip(EwaldOptions const& pmeOpts, int debugIn) {
  if (recipParams_.InitRecip(pmeOpts, debugIn)) return 1;
  PrintRecipOpts();
  return 0;
}

/** Print options to stdout. */
void PME_Recip::PrintRecipOpts() const {
  if (distKernelExponent_ == 1)
    mprintf("\tParticle Mesh Coulomb Reciprocal opts:\n");
  else if (distKernelExponent_ == 6)
    mprintf("\tParticle Mesh LJ Reciprocal opts:\n");
  else
    mprintf("\tReciprocal opts (distance kernel exponent= %i)\n", distKernelExponent_);
  recipParams_.PrintRecipOpts();
}

/** \return 1 if the box type is invalid for helPME */
int PME_Recip::set_lattice(PMEInstanceD::LatticeType& lattice, Box const& boxIn) {
  lattice = PMEInstanceD::LatticeType::XAligned;
  // TODO just pass in Ucell when helPME supports it
  //boxIn.PrintDebug("pme");
  if (!boxIn.Is_X_Aligned()) {
    if (boxIn.Is_Symmetric())
      lattice = PMEInstanceD::LatticeType::ShapeMatrix;
    else {
      mprinterr("Error: Unit cell is not X-aligned or symmetric; cannot set PME recip grid.\n");
      return 1;
    }
  }
  return 0;
}

/** \return Reciprocal space part of PME energy calc. */
// TODO currently helPME needs the coords/charge arrays to be non-const, need to fix that
double PME_Recip::Recip_ParticleMesh(Darray& coordsDin, Box const& boxIn,
                                     Darray& ChargeIn, double ew_coeffIn)
{
  t_recip_.Start();
  // This essentially makes coordsD and chargesD point to arrays.
  Mat coordsD(&coordsDin[0], ChargeIn.size(), 3);
  Mat chargesD(&ChargeIn[0], ChargeIn.size(), 1);
  int nfft1, nfft2, nfft3;

  if ( recipParams_.DetermineNfft(nfft1, nfft2, nfft3, boxIn) ) {
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
  pme_object_.setup(distKernelExponent_, ew_coeffIn, recipParams_.Order(), nfft1, nfft2, nfft3, scaleFac_, 0);
  // Check the unit cell vectors
  PMEInstanceD::LatticeType lattice;
  if (set_lattice(lattice, boxIn)) return 0;
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
                                lattice);
  //t_calc_.Start();
  double erecip = pme_object_.computeERec(0, chargesD, coordsD);
  //t_calc_.Stop();
  t_recip_.Stop();
  return erecip;
}

/** \return Reciprocal space part of PME energy calc. */
// TODO currently helPME needs the coords/charge arrays to be non-const, need to fix that
double PME_Recip::Recip_Decomp(Darray& atom_recip,
                               Darray& coordsDin, Box const& boxIn,
                               Darray& ChargeIn, double ew_coeffIn)
{
  t_recip_.Start();
  atom_recip.resize( ChargeIn.size() );
  // This essentially makes coordsD and chargesD point to arrays.
  Mat coordsD(&coordsDin[0], ChargeIn.size(), 3);
  Mat chargesD(&ChargeIn[0], ChargeIn.size(), 1);
  int nfft1, nfft2, nfft3;

  if ( recipParams_.DetermineNfft(nfft1, nfft2, nfft3, boxIn) ) {
    mprinterr("Error: Could not determine FFT grid spacing.\n");
    return 0.0;
  }

  // Instantiate double precision PME object
  pme_object_.setup(distKernelExponent_, ew_coeffIn, recipParams_.Order(), nfft1, nfft2, nfft3, scaleFac_, 0);
  // Check the unit cell vectors
  PMEInstanceD::LatticeType lattice;
  if (set_lattice(lattice, boxIn)) return 0;
  // Sets the unit cell lattice vectors, with units consistent with those used to specify coordinates.
  pme_object_.setLatticeVectors(boxIn.Param(Box::X), boxIn.Param(Box::Y), boxIn.Param(Box::Z),
                                boxIn.Param(Box::ALPHA), boxIn.Param(Box::BETA), boxIn.Param(Box::GAMMA),
                                lattice);
  //t_calc_.Start();
  // TODO precalc
  Mat e_potentialD_(ChargeIn.size(), 4);
  e_potentialD_.setConstant(0.0);
  pme_object_.computePRec(0, chargesD, coordsD, coordsD, 1, e_potentialD_);
  double erecip = 0;
  for(unsigned int i = 0; i < ChargeIn.size(); i++)
  {
    atom_recip[i]=0.5 * ChargeIn[i] * e_potentialD_(i,0);
    erecip += atom_recip[i];
  }

  //t_calc_.Stop();
  t_recip_.Stop();
  return erecip;
}
#endif /* LIBPME */

#ifdef LIBPME
#include "Ewald_ParticleMesh.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
Ewald_ParticleMesh::Ewald_ParticleMesh() : order_(6)
{
  nfft_[0] = -1;
  nfft_[1] = -1;
  nfft_[2] = -1;
}

/** \return true if given number is a product of powers of 2, 3, or 5. */
static inline bool check_prime_factors(int nIn) {
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
int Ewald_ParticleMesh::ComputeNFFT(double len) {
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
int Ewald_ParticleMesh::DetermineNfft(int& nfft1, int& nfft2, int& nfft3, Box const& boxIn) const
{
   if (nfft1 < 1) {
    // Need even dimension for X direction
    nfft1 = ComputeNFFT( (boxIn.BoxX() + 1.0) * 0.5 );
    nfft1 *= 2;
  }
  if (nfft2 < 1)
    nfft2 = ComputeNFFT( boxIn.BoxY() );
  if (nfft3 < 1)
    nfft3 = ComputeNFFT( boxIn.BoxZ() );

  if (nfft1 < 1 || nfft2 < 1 || nfft3 < 1) {
    mprinterr("Error: Bad NFFT values: %i %i %i\n", nfft1, nfft2, nfft3);
    return 1;
  }
  if (debug_ > 0) mprintf("DEBUG: NFFTs: %i %i %i\n", nfft1, nfft2, nfft3);

  return 0;
}

/** Set up PME parameters. */
int Ewald_ParticleMesh::Init(Box const& boxIn, double cutoffIn, double dsumTolIn,
                    double ew_coeffIn, double skinnbIn, double erfcTableDxIn, 
                    int orderIn, int debugIn, const int* nfftIn)
{
  if (CheckInput(boxIn, debugIn, cutoffIn, dsumTolIn, ew_coeffIn, erfcTableDxIn, skinnbIn))
    return 1;
  if (nfftIn != 0)
    std::copy(nfftIn, nfftIn+3, nfft_);
  order_ = orderIn;

  // Set defaults if necessary
  if (order_ < 1) order_ = 6;

  mprintf("\tParticle Mesh Ewald params:\n");
  mprintf("\t  Cutoff= %g   Direct Sum Tol= %g   Ewald coeff.= %g  NB skin= %g\n",
          cutoff_, dsumTol_, ew_coeff_, skinnbIn);
  mprintf("\t  Bspline order= %i\n", order_);
  mprintf("\t  Erfc table dx= %g, size= %zu\n", erfcTableDx_, erfc_table_.size()/4);
  mprintf("\t ");
  for (int i = 0; i != 3; i++)
    if (nfft_[i] == -1)
      mprintf(" NFFT%i=auto", i+1);
    else
      mprintf(" NFFT%i=%i", i+1, nfft_[i]);
  mprintf("\n");

  // Set up pair list
  Matrix_3x3 ucell, recip;
  boxIn.ToRecip(ucell, recip);
  Vec3 recipLengths = boxIn.RecipLengths(recip);
  if (Setup_Pairlist(boxIn, recipLengths, skinnbIn)) return 1;

  return 0;
}

/** Setup PME calculation. */
int Ewald_ParticleMesh::Setup(Topology const& topIn, AtomMask const& maskIn) {
  CalculateCharges(topIn, maskIn);
  coordsD_  = libpme::Mat<double>(maskIn.Nselected(), 3);
  // This essentially makes chargesD_ point to the Charge_ array.
  chargesD_ = libpme::mapMat<double>(&Charge_[0], maskIn.Nselected(), 1);
  SetupExcluded(topIn);
  return 0;
}

/** Calculate full nonbonded energy with PME */
double Ewald_ParticleMesh::CalcEnergy(Frame const& frameIn, AtomMask const& maskIn)
{
  t_total_.Start();
  Matrix_3x3 ucell, recip;
  double volume = frameIn.BoxCrd().ToRecip(ucell, recip);
  double e_self = Self( volume );

  pairList_.CreatePairList(frameIn, ucell, recip, maskIn);

  // TODO make more efficient
  int idx = 0;
  for (AtomMask::const_iterator atm = maskIn.begin(); atm != maskIn.end(); ++atm, ++idx) {
    const double* XYZ = frameIn.XYZ( *atm );
    coordsD_(idx, 0) = XYZ[0];
    coordsD_(idx, 1) = XYZ[1];
    coordsD_(idx, 2) = XYZ[2];
  }

//  MapCoords(frameIn, ucell, recip, maskIn);
  double e_recip = Recip_ParticleMesh( coordsD_, chargesD_, frameIn.BoxCrd() );
  double e_adjust = 0.0;
  double e_direct = Direct( pairList_, e_adjust );
  if (debug_ > 0)
    mprintf("DEBUG: Eself= %20.10f   Erecip= %20.10f   Edirect= %20.10f  Eadjust= %20.10f\n",
            e_self, e_recip, e_direct, e_adjust);
  t_total_.Stop();
  return e_self + e_recip + e_direct + e_adjust;
}

// Ewald::Recip_ParticleMesh()
double Ewald_ParticleMesh::Recip_ParticleMesh(libpme::Mat<double> const& coordsD,
                                 libpme::Mat<double> const& chargesD, Box const& boxIn)
{
  t_recip_.Start();
  int nfft1 = nfft_[0];
  int nfft2 = nfft_[1];
  int nfft3 = nfft_[2];
  if ( DetermineNfft(nfft1, nfft2, nfft3, boxIn) ) {
    mprinterr("Error: Could not determine grid spacing.\n");
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
  //       8 = number of nodes used for the rec space PME calculation.
  //       9 = max # threads to use for each MPI instance; 0 = all available threads used.
  // NOTE: Scale factor for Charmm is 332.0716
  // NOTE: The electrostatic constant has been baked into the Charge_ array already.
  auto pme_object = std::unique_ptr<PMEInstanceD>(new PMEInstanceD(1, ew_coeff_, order_, nfft1, nfft2, nfft3, 1.0, 1, 0));
  // Sets the unit cell lattice vectors, with units consistent with those used to specify coordinates.
  // Args: 1 = the A lattice parameter in units consistent with the coordinates.
  //       2 = the B lattice parameter in units consistent with the coordinates.
  //       3 = the C lattice parameter in units consistent with the coordinates.
  //       4 = the alpha lattice parameter in degrees.
  //       5 = the beta lattice parameter in degrees.
  //       6 = the gamma lattice parameter in degrees.
  pme_object->setLatticeVectors(boxIn.BoxX(), boxIn.BoxY(), boxIn.BoxZ(), boxIn.Alpha(), boxIn.Beta(), boxIn.Gamma());
  double erecip = pme_object->computeERec(0, chargesD, coordsD);
  t_recip_.Stop();
  return erecip;
}


#endif /* LIBPME */

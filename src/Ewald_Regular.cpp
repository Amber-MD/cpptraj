#include "Ewald_Regular.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "StringRoutines.h" // ByteString
#ifdef _OPENMP
# include <omp.h>
#endif

/// CONSTRUCTOR
Ewald_Regular::Ewald_Regular() :
# ifdef _OPENMP
  multCut_(0),
# endif
  maxexp_(0.0),
  rsumTol_(0.0),
  maxmlim_(0)
{
  mlimit_[0] = 0;
  mlimit_[1] = 0;
  mlimit_[2] = 0;
}

/** \return maxexp value based on mlimits */
double Ewald_Regular::FindMaxexpFromMlim(const int* mlimit, Matrix_3x3 const& recip) {
  double maxexp = DABS( (double)mlimit[0] * recip[0] );
  double z2     = DABS( (double)mlimit[1] * recip[4] );
  maxexp = std::max(maxexp, z2);
  double z3     = DABS( (double)mlimit[2] * recip[8] );
  maxexp = std::max(maxexp, z3);
  return maxexp;
}

/** \return maxexp value based on Ewald coefficient and reciprocal sum tolerance. */
double Ewald_Regular::FindMaxexpFromTol(double ewCoeff, double rsumTol) {
  double xval = 0.5;
  int nloop = 0;
  double term = 0.0;
  do {
    xval = 2.0 * xval;
    nloop++;
    double yval = Constants::PI * xval / ewCoeff;
    term = 2.0 * ewCoeff * erfc_func(yval) * INVSQRTPI_;
  } while (term >= rsumTol);

  // Binary search tolerance is 2^-60
  int ntimes = nloop + 60;
  double xlo = 0.0;
  double xhi = xval;
  for (int i = 0; i != ntimes; i++) {
    xval = (xlo + xhi) / 2.0;
    double yval = Constants::PI * xval / ewCoeff;
    double term = 2.0 * ewCoeff * erfc_func(yval) * INVSQRTPI_;
    if (term > rsumTol)
      xlo = xval;
    else
      xhi = xval;
  }
  mprintf("\tMaxExp for Ewald coefficient %g, direct sum tol %g is %g\n",
          ewCoeff, rsumTol, xval);
  return xval;
}

static inline int IABS(int    xIn) { if (xIn < 0  ) return -xIn; else return xIn; }

/** Get mlimits. */
void Ewald_Regular::GetMlimits(int* mlimit, double maxexp, double eigmin, 
                       Vec3 const& reclng, Matrix_3x3 const& recip)
{
  //mprintf("DEBUG: Recip lengths %12.4f%12.4f%12.4f\n", reclng[0], reclng[1], reclng[2]);

  int mtop1 = (int)(reclng[0] * maxexp / sqrt(eigmin));
  int mtop2 = (int)(reclng[1] * maxexp / sqrt(eigmin));
  int mtop3 = (int)(reclng[2] * maxexp / sqrt(eigmin));

  int nrecvecs = 0;
  mlimit[0] = 0;
  mlimit[1] = 0;
  mlimit[2] = 0;
  double maxexp2 = maxexp * maxexp;
  for (int m1 = -mtop1; m1 <= mtop1; m1++) {
    for (int m2 = -mtop2; m2 <= mtop2; m2++) {
      for (int m3 = -mtop3; m3 <= mtop3; m3++) {
        Vec3 Zvec = recip.TransposeMult( Vec3(m1,m2,m3) );
        if ( Zvec.Magnitude2() <= maxexp2 ) {
          nrecvecs++;
          mlimit[0] = std::max( mlimit[0], IABS(m1) );
          mlimit[1] = std::max( mlimit[1], IABS(m2) );
          mlimit[2] = std::max( mlimit[2], IABS(m3) );
        }
      }
    }
  }
  mprintf("\tNumber of reciprocal vectors: %i\n", nrecvecs);
}

/** Init regular Ewald calculation. */
int Ewald_Regular::Init(Box const& boxIn, double cutoffIn, double dsumTolIn, double rsumTolIn,
                     double ew_coeffIn, double maxexpIn, double skinnbIn,
                     double erfcTableDxIn, int debugIn, const int* mlimitsIn)
{
  if (CheckInput(boxIn, debugIn, cutoffIn, dsumTolIn, ew_coeffIn, erfcTableDxIn, skinnbIn))
    return 1;
  rsumTol_ = rsumTolIn;
  maxexp_ = maxexpIn;
  Matrix_3x3 ucell, recip;
  boxIn.ToRecip(ucell, recip);
  if (mlimitsIn != 0)
    std::copy(mlimitsIn, mlimitsIn+3, mlimit_);

  // Check input
  if (mlimit_[0] < 0 || mlimit_[1] < 0 || mlimit_[2] < 0) {
    mprinterr("Error: Cannot specify negative mlimit values.\n");
    return 1;
  }
  maxmlim_ = mlimit_[0];
  maxmlim_ = std::max(maxmlim_, mlimit_[1]);
  maxmlim_ = std::max(maxmlim_, mlimit_[2]);
  if (maxexp_ < 0.0) {
    mprinterr("Error: maxexp is less than 0.0\n");
    return 1;
  }

  // Set defaults if necessary
  if (rsumTol_ < Constants::SMALL)
    rsumTol_ = 5E-5;
  Vec3 recipLengths = boxIn.RecipLengths(recip);
  if (maxmlim_ > 0)
    maxexp_ = FindMaxexpFromMlim(mlimit_, recip);
  else {
    if ( maxexp_ < Constants::SMALL )
      maxexp_ = FindMaxexpFromTol(ew_coeff_, rsumTol_);
    // eigmin typically bigger than this unless cell is badly distorted.
    double eigmin = 0.5;
    // Calculate lengths of reciprocal vectors
    GetMlimits(mlimit_, maxexp_, eigmin, recipLengths, recip);
    maxmlim_ = mlimit_[0];
    maxmlim_ = std::max(maxmlim_, mlimit_[1]);
    maxmlim_ = std::max(maxmlim_, mlimit_[2]);
  }

  mprintf("\tEwald params:\n");
  mprintf("\t  Cutoff= %g   Direct Sum Tol= %g   Ewald coeff.= %g\n",
          cutoff_, dsumTol_, ew_coeff_);
  mprintf("\t  MaxExp= %g   Recip. Sum Tol= %g   NB skin= %g\n",
          maxexp_, rsumTol_, skinnbIn);
  mprintf("\t  Erfc table dx= %g, size= %zu\n", erfcTableDx_, erfc_table_.size()/4);
  mprintf("\t  mlimits= {%i,%i,%i} Max=%i\n", mlimit_[0], mlimit_[1], mlimit_[2], maxmlim_);
  // Set up pair list
  if (Setup_Pairlist(boxIn, recipLengths, skinnbIn)) return 1;

  return 0;
}

/** Setup regular Ewald calculation. */
int Ewald_Regular::Setup(Topology const& topIn, AtomMask const& maskIn) {
  CalculateCharges(topIn, maskIn);

  // Build exponential factors for use in structure factors.
  // These arrays are laid out in 1D; value for each atom at each m, i.e.
  // A0M0 A1M0 A2M0 ... ANM0 A0M1 ... ANMX
  // Number of M values is the max + 1.
  int mmax = maxmlim_ + 1;
  unsigned int tsize = maskIn.Nselected() * mmax;
  cosf1_.assign( tsize, 1.0 );
  cosf2_.assign( tsize, 1.0 );
  cosf3_.assign( tsize, 1.0 );
  sinf1_.assign( tsize, 0.0 );
  sinf2_.assign( tsize, 0.0 );
  sinf3_.assign( tsize, 0.0 );
  mprintf("\tMemory used by trig tables: %s\n",
          ByteString(6*tsize*sizeof(double), BYTE_DECIMAL).c_str());
  // M0
//  for (int i = 0; i != maskIn.Nselected(); i++) {
//    cosf1_.push_back( 1.0 );
//    cosf2_.push_back( 1.0 );
//    cosf3_.push_back( 1.0 );
//    sinf1_.push_back( 0.0 );
//    sinf2_.push_back( 0.0 );
//    sinf3_.push_back( 0.0 );
// }

  SetupExcluded(topIn);

# ifdef _OPENMP
  // Pre-calculate m1 and m2 indices
  mlim1_.clear();
  mlim2_.clear();
  multCut_ = 0;
  for (int m1 = 0; m1 <= mlimit_[0]; m1++) {
    for (int m2 = -mlimit_[1]; m2 <= mlimit_[1]; m2++) {
      mlim1_.push_back( m1 );
      mlim2_.push_back( m2 );
    }
    // After this index (end of m1 == 0) multiplier must be 2.0
    if (m1 == 0)
      multCut_ = (int)mlim1_.size();
  }
  // Each thread will need its own space for trig math
  int numthreads;
# pragma omp parallel
  {
#   pragma omp master
    {
      numthreads = omp_get_num_threads();
      mprintf("\tParallelizing calculation with %i threads\n", numthreads);
    }
  }
  unsigned int asize = (unsigned int)maskIn.Nselected() * (unsigned int)numthreads;
  c12_.resize( asize );
  s12_.resize( asize );
  c3_.resize(  asize );
  s3_.resize(  asize );
#else
  c12_.resize( maskIn.Nselected() );
  s12_.resize( maskIn.Nselected() );
  c3_.resize(  maskIn.Nselected() );
  s3_.resize(  maskIn.Nselected() );
# endif
  return 0;
}

/** Reciprocal space energy counteracting the neutralizing charge distribution. */
double Ewald_Regular::Recip_Regular(Matrix_3x3 const& recip, double volume) {
  t_recip_.Start();
  double fac = (Constants::PI*Constants::PI) / (ew_coeff_ * ew_coeff_);
  double maxexp2 = maxexp_ * maxexp_;
  double ene = 0.0;
  Varray const& Frac = pairList_.FracCoords();
  // Number of M values is the max + 1.
  int mmax = maxmlim_ + 1;
  // Build exponential factors for use in structure factors.
  // These arrays are laid out in 1D; value for each atom at each m, i.e.
  // A0M0 A1M0 A2M0 ... ANM0 A0M1 ... ANMX
  // M0 is done in EwaldSetup()
  t_trig_tables_.Start();
  unsigned int mnidx = Frac.size();
  // M1
  for (unsigned int i = 0; i != Frac.size(); i++, mnidx++) {
    //mprintf("FRAC: %6i%20.10f%20.10f%20.10f\n", i+1, Frac[i][0], Frac[i][1], Frac[i][2]);
    cosf1_[mnidx] = cos(Constants::TWOPI * Frac[i][0]);
    cosf2_[mnidx] = cos(Constants::TWOPI * Frac[i][1]);
    cosf3_[mnidx] = cos(Constants::TWOPI * Frac[i][2]);
    sinf1_[mnidx] = sin(Constants::TWOPI * Frac[i][0]);
    sinf2_[mnidx] = sin(Constants::TWOPI * Frac[i][1]);
    sinf3_[mnidx] = sin(Constants::TWOPI * Frac[i][2]);
  }
  // M2-MX
  // Get the higher factors by recursion using trig addition rules.
  // Negative values of M by complex conjugation, or even cosf, odd sinf.
  // idx will always point to M-1 values
  unsigned int idx = Frac.size();
  for (int m = 2; m < mmax; m++) {
    // Set m1idx to beginning of M1 values.
    unsigned int m1idx = Frac.size();
    for (unsigned int i = 0; i != Frac.size(); i++, idx++, m1idx++, mnidx++) {
      cosf1_[mnidx] = cosf1_[idx]*cosf1_[m1idx] - sinf1_[idx]*sinf1_[m1idx];
      cosf2_[mnidx] = cosf2_[idx]*cosf2_[m1idx] - sinf2_[idx]*sinf2_[m1idx];
      cosf3_[mnidx] = cosf3_[idx]*cosf3_[m1idx] - sinf3_[idx]*sinf3_[m1idx];
      sinf1_[mnidx] = sinf1_[idx]*cosf1_[m1idx] + cosf1_[idx]*sinf1_[m1idx];
      sinf2_[mnidx] = sinf2_[idx]*cosf2_[m1idx] + cosf2_[idx]*sinf2_[m1idx];
      sinf3_[mnidx] = sinf3_[idx]*cosf3_[m1idx] + cosf3_[idx]*sinf3_[m1idx];
    }
  }
  // DEBUG
/*  unsigned int midx = 0;
  for (int m = 0; m != mmax; m++) {
    for (unsigned int i = 0; i != Frac.size(); i++, midx++)
      mprintf("TRIG: %6i%6u%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n", m,i+1,
               cosf1_[midx], cosf2_[midx], cosf3_[midx],
               sinf1_[midx], sinf2_[midx], sinf3_[midx]);
  }*/
  t_trig_tables_.Stop();
# ifdef _OPENMP
  double mult;
  unsigned int offset;
  int mlim_idx;
  int mlim_end = (int)mlim1_.size();
  double *c12, *s12, *c3, *s3;
# pragma omp parallel private(mult,mlim_idx,c12,s12,c3,s3,offset) reduction(+:ene)
  {
  offset = (unsigned int)omp_get_thread_num() * Frac.size();
  c12 = &c12_[0] + offset;
  s12 = &s12_[0] + offset;
  c3  = &c3_[0]  + offset;
  s3  = &s3_[0]  + offset;
# pragma omp for
  for (mlim_idx = 0; mlim_idx < mlim_end; mlim_idx++)
  {
      if (mlim_idx < multCut_)
        mult = 1.0;
      else
        mult = 2.0;
      int m1 = mlim1_[mlim_idx];
      int m2 = mlim2_[mlim_idx];
# else
  Darray& c12 = c12_;
  Darray& s12 = s12_;
  Darray& c3 = c3_;
  Darray& s3 = s3_;
  double mult = 1.0;
  for (int m1 = 0; m1 <= mlimit_[0]; m1++)
  {
    for (int m2 = -mlimit_[1]; m2 <= mlimit_[1]; m2++)
    {
# endif
      int m1idx = Frac.size() * m1;
      int m2idx = Frac.size() * IABS(m2);

      if (m2 < 0) {
        for (unsigned int i = 0; i != Frac.size(); i++, m1idx++, m2idx++) {
          c12[i] = cosf1_[m1idx]*cosf2_[m2idx] + sinf1_[m1idx]*sinf2_[m2idx];
          s12[i] = sinf1_[m1idx]*cosf2_[m2idx] - cosf1_[m1idx]*sinf2_[m2idx];
        }
      } else {
        for (unsigned int i = 0; i != Frac.size(); i++, m1idx++, m2idx++) {
          c12[i] = cosf1_[m1idx]*cosf2_[m2idx] - sinf1_[m1idx]*sinf2_[m2idx];
          s12[i] = sinf1_[m1idx]*cosf2_[m2idx] + cosf1_[m1idx]*sinf2_[m2idx];
        }
      }
      for (int m3 = -mlimit_[2]; m3 <= mlimit_[2]; m3++)
      {
        // Columns of recip are reciprocal unit cell vecs, so
        // mhat contains Cartesian components of recip vector M.
        Vec3 mhat = recip.TransposeMult( Vec3(m1, m2, m3) );
        double msq = mhat.Magnitude2();
        double denom = Constants::PI * volume * msq;
        double eterm = 0.0;
//        double vterm = 0.0;
        if ( m1*m1 + m2*m2 + m3*m3 > 0 ) {
          eterm = exp(-fac*msq) / denom;
//          vterm = 2.0 * (fac*msq + 1.0) / msq;
        }
        // mult takes care to double count for symmetry. Can take care of
        // with eterm.
        eterm *= mult;
        if (msq < maxexp2) {
          int m3idx = Frac.size() * IABS(m3);
          // Get the product of complex exponentials.
          if (m3 < 0) {
            for (unsigned int i = 0; i != Frac.size(); i++, m3idx++) {
              c3[i] = c12[i]*cosf3_[m3idx] + s12[i]*sinf3_[m3idx];
              s3[i] = s12[i]*cosf3_[m3idx] - c12[i]*sinf3_[m3idx];
            }
          } else {
            for (unsigned int i = 0; i != Frac.size(); i++, m3idx++) {
              c3[i] = c12[i]*cosf3_[m3idx] - s12[i]*sinf3_[m3idx];
              s3[i] = s12[i]*cosf3_[m3idx] + c12[i]*sinf3_[m3idx];
            }
          }
          // Get the structure factor
          double cstruct = 0.0;
          double sstruct = 0.0;
          for (unsigned int i = 0; i != Frac.size(); i++) {
            cstruct += Charge_[i] * c3[i];
            sstruct += Charge_[i] * s3[i];
          }
          double struc2 = cstruct*cstruct + sstruct*sstruct;
          ene += eterm * struc2;
          //mprintf("LOOP: %3i%3i%3i ENE= %20.10f\n", m1, m2, m3, ene);
        } // END IF msq < maxexp2
      } // END loop over m3
# ifdef _OPENMP
  } // END loop over mlim_idx
  } // END pragma omp parallel
# else
    } // END loop over m2
    mult = 2.0;
  } // END loop over m1
# endif
  t_recip_.Stop();
  return ene * 0.5;
}

/** Calculate Ewald energy. Faster version that uses pair list. */
double Ewald_Regular::CalcEnergy(Frame const& frameIn, AtomMask const& maskIn, double& e_vdw)
{
  t_total_.Start();
  Matrix_3x3 ucell, recip;
  double volume = frameIn.BoxCrd().ToRecip(ucell, recip);
  double e_self = Self( volume );
  double e_vdwr = Vdw_Correction( volume );

  pairList_.CreatePairList(frameIn, ucell, recip, maskIn);

//  MapCoords(frameIn, ucell, recip, maskIn);
  double e_recip = Recip_Regular( recip, volume );
  double e_adjust = 0.0;
         e_vdw = 0.0;
  double e_direct = Direct( pairList_, e_adjust, e_vdw );
  if (debug_ > 0)
    mprintf("DEBUG: Eself= %20.10f   Erecip= %20.10f   Edirect= %20.10f  Eadjust= %20.10f  Evdw= %20.10f\n",
            e_self, e_recip, e_direct, e_adjust, e_vdw);
  e_vdw += e_vdwr;
  t_total_.Stop();
  return e_self + e_recip + e_direct + e_adjust;
}


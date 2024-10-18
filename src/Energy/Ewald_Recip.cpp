#include "Ewald_Recip.h"
#include "ErfcFxn.h" // erfc_func
#include "EwaldParams.h" // INVSQRTPI
#include "../Box.h"
#include "../Constants.h"
#include "../CpptrajStdio.h"
#include "../EwaldOptions.h"
#include "../Matrix_3x3.h"
#include "../StringRoutines.h" // ByteString
#include "../Vec3.h"
#include <cmath> // sqrt, cos, sin
#include <algorithm> // std::max
#ifdef _OPENMP
# include <omp.h>
#endif 

using namespace Cpptraj::Energy;

// Absolute value functions
static inline int IABS(   int    xIn) { if (xIn < 0  ) return -xIn; else return xIn; }
static inline double DABS(double xIn) { if (xIn < 0  ) return -xIn; else return xIn; }

/** CONSTRUCTOR */
Ewald_Recip::Ewald_Recip() :
# ifdef _OPENMP
  multCut_(0),
# endif
  fac_(0.0),
  maxexp_(0.0),
  rsumTol_(0.0),
  maxmlim_(0)
{
  mlimit_[0] = 0;
  mlimit_[1] = 0;
  mlimit_[2] = 0;
}

/** \return maxexp value based on mlimits */
double Ewald_Recip::FindMaxexpFromMlim(const int* mlimit, Matrix_3x3 const& recip) {
  double maxexp = DABS( (double)mlimit[0] * recip[0] );
  double z2     = DABS( (double)mlimit[1] * recip[4] );
  maxexp = std::max(maxexp, z2);
  double z3     = DABS( (double)mlimit[2] * recip[8] );
  maxexp = std::max(maxexp, z3);
  return maxexp;
}

/** \return maxexp value based on Ewald coefficient and reciprocal sum tolerance. */
double Ewald_Recip::FindMaxexpFromTol(double ewCoeff, double rsumTol) {
  double xval = 0.5;
  int nloop = 0;
  double term = 0.0;
  do {
    xval = 2.0 * xval;
    nloop++;
    double yval = Constants::PI * xval / ewCoeff;
    term = 2.0 * ewCoeff * ErfcFxn::erfc_func(yval) * EwaldParams::INVSQRTPI();
  } while (term >= rsumTol);

  // Binary search tolerance is 2^-60
  int ntimes = nloop + 60;
  double xlo = 0.0;
  double xhi = xval;
  for (int i = 0; i != ntimes; i++) {
    xval = (xlo + xhi) / 2.0;
    double yval = Constants::PI * xval / ewCoeff;
    double term = 2.0 * ewCoeff * ErfcFxn::erfc_func(yval) * EwaldParams::INVSQRTPI();
    if (term > rsumTol)
      xlo = xval;
    else
      xhi = xval;
  }
  mprintf("\t  MaxExp determined from Ewald coefficient %g, direct sum tol %g is %g\n",
          ewCoeff, rsumTol, xval);
  return xval;
}


/** Get mlimits. */
void Ewald_Recip::GetMlimits(int* mlimit, double maxexp, double eigmin, 
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
  mprintf("\t  Number of reciprocal vectors: %i\n", nrecvecs);
}

/** Init */
int Ewald_Recip::InitRecip(EwaldOptions const& ewOpts, double ew_coeffIn,
                           Box const& boxIn, int debugIn)
{
  debug_ = debugIn;
  rsumTol_ = ewOpts.RsumTol();
  maxexp_ = ewOpts.MaxExp();
  fac_ = (Constants::PI*Constants::PI) / (ew_coeffIn * ew_coeffIn);
  mlimit_[0] = ewOpts.Mlimits1();
  mlimit_[1] = ewOpts.Mlimits2();
  mlimit_[2] = ewOpts.Mlimits3();

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
  mprintf("\tRecip opts (regular Ewald):\n");
  if (rsumTol_ < Constants::SMALL)
    rsumTol_ = 5E-5;
  if (maxmlim_ > 0)
    maxexp_ = FindMaxexpFromMlim(mlimit_, boxIn.FracCell());
  else {
    if ( maxexp_ < Constants::SMALL )
      maxexp_ = FindMaxexpFromTol(ew_coeffIn, rsumTol_);
    // eigmin typically bigger than this unless cell is badly distorted.
    double eigmin = 0.5;
    // Calculate lengths of reciprocal vectors
    GetMlimits(mlimit_, maxexp_, eigmin, boxIn.RecipLengths(), boxIn.FracCell());
    maxmlim_ = mlimit_[0];
    maxmlim_ = std::max(maxmlim_, mlimit_[1]);
    maxmlim_ = std::max(maxmlim_, mlimit_[2]);
  }

  PrintRecipOpts();
  return 0;
}

/** Print options to stdout */
void Ewald_Recip::PrintRecipOpts() const {
  mprintf("\t  MaxExp= %g   Recip. Sum Tol= %g\n", maxexp_, rsumTol_);
  //mprintf("\t  Erfc table dx= %g, size= %zu\n", erfcTableDx_, erfc_table_.size()/4);
  mprintf("\t  mlimits= {%i,%i,%i} Max=%i\n", mlimit_[0], mlimit_[1], mlimit_[2], maxmlim_);
}

/** Setup trig tables for given number of selected atoms. */
int Ewald_Recip::SetupRecip(int nselected) {
  // Build exponential factors for use in structure factors.
  // These arrays are laid out in 1D; value for each atom at each m, i.e.
  // A0M0 A1M0 A2M0 ... ANM0 A0M1 ... ANMX
  // Number of M values is the max + 1.
  int mmax = maxmlim_ + 1;
  unsigned int tsize = nselected * mmax;
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
  unsigned int asize = (unsigned int)nselected * (unsigned int)numthreads;
  c12_.resize( asize );
  s12_.resize( asize );
  c3_.resize(  asize );
  s3_.resize(  asize );
#else /* _OPENMP */
  c12_.resize( nselected );
  s12_.resize( nselected );
  c3_.resize(  nselected );
  s3_.resize(  nselected );
# endif /* _OPENMP */

  return 0;
}

/** Reciprocal space energy counteracting the neutralizing charge distribution. */
double Ewald_Recip::Recip_Regular(Matrix_3x3 const& recip, double volume,
                                  Varray const& Frac, Darray const& Charge)
{
  t_recip_.Start();
  double maxexp2 = maxexp_ * maxexp_;
  double ene = 0.0;
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
          eterm = exp(-fac_*msq) / denom;
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
            cstruct += Charge[i] * c3[i];
            sstruct += Charge[i] * s3[i];
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

void Ewald_Recip::PrintTiming(double total) const {
  t_recip_.WriteTiming(2,  "Recip:     ", total);
  if (t_trig_tables_.Total() > 0.0)
    t_trig_tables_.WriteTiming(3, "Calc trig tables:", t_recip_.Total());
}

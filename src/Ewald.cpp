#include <cmath>
#include <algorithm> // std::min, std::max
#include "Ewald.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "StringRoutines.h"
#include "Spline.h"
#ifdef _OPENMP
# include <omp.h>
#endif

Ewald::Ewald() :
  sumq_(0.0),
  sumq2_(0.0),
  ew_coeff_(0.0),
  maxexp_(0.0),
  cutoff_(0.0),
  dsumTol_(0.0),
  rsumTol_(0.0),
  erfcTableDx_(0.0),
  one_over_Dx_(0.0),
  maxmlim_(0),
  debug_(0)
{
  mlimit_[0] = 0;
  mlimit_[1] = 0;
  mlimit_[2] = 0;
# ifdef DEBUG_EWALD
  // Save fractional translations for 1 cell in each direction (and primary cell).
  // This is only for the non-pairlist version of direct.
  Cells_.reserve( 27 );
  for (int ix = -1; ix < 2; ix++)
    for (int iy = -1; iy < 2; iy++)
      for (int iz = -1; iz < 2; iz++)
        Cells_.push_back( Vec3(ix, iy, iz) );
# endif
}

const double Ewald::INVSQRTPI_ = 1.0 / sqrt(Constants::PI);

static inline double DABS(double xIn) { if (xIn < 0.0) return -xIn; else return xIn; }
static inline int    IABS(int    xIn) { if (xIn < 0  ) return -xIn; else return xIn; }

/** Complimentary error function: 2/sqrt(PI) * SUM[exp(-t^2)*dt]
  * Original code: SANDER: erfcfun.F90
  */
double Ewald::erfc_func(double xIn) {
  double erfc;
  double absx = DABS( xIn );
    
  if (xIn > 26.0)
    erfc = 0.0;
  else if (xIn < -5.5)
    erfc = 2.0;
  else if (absx <= 0.5) {
    double cval = xIn * xIn;
    double pval = ((-0.356098437018154E-1*cval+0.699638348861914E1)*cval + 0.219792616182942E2) *
                  cval + 0.242667955230532E3;
    double qval = ((cval+0.150827976304078E2)*cval+0.911649054045149E2)*cval + 0.215058875869861E3;
    double erf = xIn * pval/qval;
    erfc = 1.0 - erf;
  } else if (absx < 4.0) {
    double cval = absx;
    double pval=((((((-0.136864857382717E-6*cval+0.564195517478974)*cval+
                     0.721175825088309E1)*cval+0.431622272220567E2)*cval+
                   0.152989285046940E3)*cval+0.339320816734344E3)*cval+
                 0.451918953711873E3)*cval+0.300459261020162E3;
    double qval=((((((cval+0.127827273196294E2)*cval+0.770001529352295E2)*cval+
                    0.277585444743988E3)*cval+0.638980264465631E3)*cval+
                  0.931354094850610E3)*cval+0.790950925327898E3)*cval+
                0.300459260956983E3;
    double nonexperfc;
    if ( xIn > 0.0 )
      nonexperfc = pval/qval;
    else
      nonexperfc = 2.0*exp(xIn*xIn) - pval/qval;
    erfc = exp(-absx*absx)*nonexperfc;
  } else {
    double cval = 1.0/(xIn*xIn);
    double pval = (((0.223192459734185E-1*cval+0.278661308609648)*cval+
                    0.226956593539687)*cval+0.494730910623251E-1)*cval+
                  0.299610707703542E-2;
    double qval = (((cval+0.198733201817135E1)*cval+0.105167510706793E1)*cval+
                   0.191308926107830)*cval+0.106209230528468E-1;
    cval = (-cval*pval/qval + 0.564189583547756)/absx;
    double nonexperfc;
    if ( xIn > 0.0 )
      nonexperfc = cval;
    else
      nonexperfc = 2.0*exp(xIn*xIn) - cval;
    erfc = exp(-absx*absx)*nonexperfc;
  }
  return erfc;
}

/** Determine Ewald coefficient from cutoff and direct sum tolerance.
  * Original Code: SANDER: findewaldcof
  */
double Ewald::FindEwaldCoefficient(double cutoff, double dsum_tol)
{
  // First get direct sum tolerance. How big must the Ewald coefficient be to
  // get terms outside the cutoff below tolerance?
  double xval = 0.5;
  int nloop = 0;
  double term = 0.0;
  do {
    xval = 2.0 * xval;
    nloop++;
    double yval = xval * cutoff;
    term = erfc_func(yval) / cutoff;
  } while (term >= dsum_tol);

  // Binary search tolerance is 2^-50
  int ntimes = nloop + 50;
  double xlo = 0.0;
  double xhi = xval;
  for (int i = 0; i != ntimes; i++) {
    xval = (xlo + xhi) / 2.0;
    double yval = xval * cutoff;
    double term = erfc_func(yval) / cutoff;
    if (term >= dsum_tol)
      xlo = xval;
    else
      xhi = xval;
  }
  mprintf("\tEwald coefficient for cut=%g, direct sum tol=%g is %g\n",
          cutoff, dsum_tol, xval);
  return xval;
}

/** \return maxexp value based on mlimits */
double Ewald::FindMaxexpFromMlim(const int* mlimit, Matrix_3x3 const& recip) {
  double maxexp = DABS( (double)mlimit[0] * recip[0] );
  double z2     = DABS( (double)mlimit[1] * recip[4] );
  maxexp = std::max(maxexp, z2);
  double z3     = DABS( (double)mlimit[2] * recip[8] );
  maxexp = std::max(maxexp, z3);
  return maxexp;
}

/** \return maxexp value based on Ewald coefficient and reciprocal sum tolerance. */
double Ewald::FindMaxexpFromTol(double ewCoeff, double rsumTol) {
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

/** Get mlimits. */
void Ewald::GetMlimits(int* mlimit, double maxexp, double eigmin, 
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

// Ewald::FillErfcTable()
void Ewald::FillErfcTable(double cutoffIn, double dxdr) {
  one_over_Dx_ = 1.0 / erfcTableDx_;
  unsigned int erfcTableSize = (unsigned int)(dxdr * one_over_Dx_ * cutoffIn * 1.5);
  Darray erfc_X, erfc_Y;
  erfc_X.reserve( erfcTableSize );
  erfc_Y.reserve( erfcTableSize );
  // Save X and Y values so we can calc the spline coefficients
  double xval = 0.0;
  for (unsigned int i = 0; i != erfcTableSize; i++) {
    double yval = erfc_func( xval );
    erfc_X.push_back( xval );
    erfc_Y.push_back( yval );
    xval += erfcTableDx_;
  }
  Spline cspline;
  cspline.CubicSpline_Coeff(erfc_X, erfc_Y);
  erfc_X.clear();
  // Store values in Spline table
  erfc_table_.reserve( erfcTableSize * 4 ); // Y B C D
  for (unsigned int i = 0; i != erfcTableSize; i++) {
    erfc_table_.push_back( erfc_Y[i] );
    erfc_table_.push_back( cspline.B_coeff()[i] );
    erfc_table_.push_back( cspline.C_coeff()[i] );
    erfc_table_.push_back( cspline.D_coeff()[i] );
  }
  // Memory saved Y values plus spline B, C, and D coefficient arrays.
  mprintf("\tMemory used by Erfc table and splines: %s\n",
          ByteString(erfc_table_.size() * sizeof(double), BYTE_DECIMAL).c_str());
}

// Ewald::ERFC()
double Ewald::ERFC(double xIn) const {
  int xidx = ((int)(one_over_Dx_ * xIn));
  double dx = xIn - ((double)xidx * erfcTableDx_);
  xidx *= 4;
  return erfc_table_[xidx] + 
         dx*(erfc_table_[xidx+1] + dx*(erfc_table_[xidx+2] + dx*erfc_table_[xidx+3]));
}

// -----------------------------------------------------------------------------
/** Set up parameters. */
int Ewald::EwaldInit(Box const& boxIn, double cutoffIn, double dsumTolIn, double rsumTolIn,
                     double ew_coeffIn, double maxexpIn, double skinnbIn,
                     double erfcTableDxIn, int debugIn, const int* mlimitsIn)
{
  debug_ = debugIn;
  cutoff_ = cutoffIn;
  dsumTol_ = dsumTolIn;
  rsumTol_ = rsumTolIn;
  ew_coeff_ = ew_coeffIn;
  maxexp_ = maxexpIn;
  erfcTableDx_ = erfcTableDxIn;
  Matrix_3x3 ucell, recip;
  boxIn.ToRecip(ucell, recip);
  if (mlimitsIn != 0)
    std::copy(mlimitsIn, mlimitsIn+3, mlimit_);

  // Check input
  if (cutoff_ < Constants::SMALL) {
    mprinterr("Error: Direct space cutoff (%g) is too small.\n", cutoff_);
    return 1;
  }
  char dir[3] = {'X', 'Y', 'Z'};
  for (int i = 0; i < 3; i++) {
    if (cutoff_ > boxIn[i]/2.0) {
      mprinterr("Error: Cutoff must be less than half the box length (%g > %g, %c)\n",
                cutoff_, boxIn[i]/2.0, dir[i]);
      return 1;
    }
  }
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
  if (skinnbIn < 0.0) {
    mprinterr("Error: skinnb is less than 0.0\n");
    return 1;
  }

  // Set defaults if necessary
  if (dsumTol_ < Constants::SMALL)
    dsumTol_ = 1E-5;
  if (rsumTol_ < Constants::SMALL)
    rsumTol_ = 5E-5;
  Vec3 recipLengths = boxIn.RecipLengths(recip);
  if (DABS(ew_coeff_) < Constants::SMALL)
    ew_coeff_ = FindEwaldCoefficient( cutoff_, dsumTol_ );
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
  if (erfcTableDx_ <= 0.0) erfcTableDx_ = 1.0 / 5000;
  // TODO make this optional
  FillErfcTable( cutoff_, ew_coeff_ );

  mprintf("\tEwald params:\n");
  mprintf("\t  Cutoff= %g   Direct Sum Tol= %g   Ewald coeff.= %g\n",
          cutoff_, dsumTol_, ew_coeff_);
  mprintf("\t  MaxExp= %g   Recip. Sum Tol= %g   NB skin= %g\n",
          maxexp_, rsumTol_, skinnbIn);
  mprintf("\t  Erfc table dx= %g, size= %zu\n", erfcTableDx_, erfc_table_.size()/4);
  mprintf("\t  mlimits= {%i,%i,%i} Max=%i\n", mlimit_[0], mlimit_[1], mlimit_[2], maxmlim_);
  // Set up pair list
  if (pairList_.InitPairList(cutoff_, skinnbIn, debugIn)) return 1;
  if (pairList_.SetupPairList( boxIn.Type(), recipLengths )) return 1;

  return 0;
}

/** Convert charges to Amber units. Calculate sum of charges and squared charges. */
void Ewald::EwaldSetup(Topology const& topIn, AtomMask const& maskIn) {
  sumq_ = 0.0;
  sumq2_ = 0.0;
  Charge_.clear();
  for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom) {
    double qi = topIn[*atom].Charge() * Constants::ELECTOAMBER;
    Charge_.push_back(qi);
    sumq_ += qi;
    sumq2_ += (qi * qi);
  }
  //mprintf("DEBUG: sumq= %20.10f   sumq2= %20.10f\n", sumq_, sumq2_);
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
  // Set up full exclusion lists.
  Excluded_.clear();
  Excluded_.resize( topIn.Natom() );
  for (int at = 0; at != topIn.Natom(); at++) {
    // Always exclude self
    Excluded_[at].insert( at );
    for (Atom::excluded_iterator excluded_atom = topIn[at].excludedbegin();
                                 excluded_atom != topIn[at].excludedend();
                               ++excluded_atom)
    {
      Excluded_[at            ].insert( *excluded_atom );
      Excluded_[*excluded_atom].insert( at             );
    }
  }
  unsigned int ex_size = 0;
  for (Iarray2D::const_iterator it = Excluded_.begin(); it != Excluded_.end(); ++it)
    ex_size += it->size();
  mprintf("\tMemory used by full exclusion list: %s\n",
          ByteString(ex_size * sizeof(int), BYTE_DECIMAL).c_str());
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
}

/** Self energy. This is the cancelling Gaussian plus the "neutralizing plasma". */
double Ewald::Self(double volume) {
  t_self_.Start();
  double d0 = -ew_coeff_ * INVSQRTPI_;
  double ene = sumq2_ * d0;
//  mprintf("DEBUG: d0= %20.10f   ene= %20.10f\n", d0, ene);
  double factor = Constants::PI / (ew_coeff_ * ew_coeff_ * volume);
  double ee_plasma = -0.5 * factor * sumq_ * sumq_;
  ene += ee_plasma;
  t_self_.Stop();
  return ene;
}

/** Reciprocal space energy counteracting the neutralizing charge distribution. */
double Ewald::Recip_Regular(Matrix_3x3 const& recip, double volume) {
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

// Ewald::Adjust()
# ifdef _OPENMP
double Ewald::Adjust(double q0, double q1, double rij) const {
  double erfc = ERFC(ew_coeff_ * rij);
  double d0 = (erfc - 1.0) / rij;
  return (q0 * q1 * d0);
}
# else
double Ewald::Adjust(double q0, double q1, double rij) {
  t_adjust_.Start();
  t_erfc_.Start();
  //double erfc = erfc_func(ew_coeff_ * rij);
  double erfc = ERFC(ew_coeff_ * rij);
  t_erfc_.Stop();
  double d0 = (erfc - 1.0) / rij;
  t_adjust_.Stop();
  return (q0 * q1 * d0);
}
# endif

//  Ewald::Direct()
/** Calculate direct space energy. This is the faster version that uses
  * a pair list. Also calculate the energy adjustment for excluded
  * atoms.
  */
double Ewald::Direct(PairList const& PL, double& e_adjust_out)
{
  t_direct_.Start();
  double cut2 = cutoff_ * cutoff_;
  double Eelec = 0.0;
  double e_adjust = 0.0;
  int cidx;
# ifdef _OPENMP
# pragma omp parallel private(cidx) reduction(+: Eelec, e_adjust)
  {
# pragma omp for
# endif
  for (cidx = 0; cidx < PL.NGridMax(); cidx++)
  {
    PairList::CellType const& thisCell = PL.Cell( cidx );
    if (thisCell.NatomsInGrid() > 0)
    {
      // cellList contains this cell index and all neighbors.
      PairList::Iarray const& cellList = thisCell.CellList();
      // transList contains index to translation for the neighbor.
      PairList::Iarray const& transList = thisCell.TransList();
      // Loop over all atoms of thisCell.
      for (PairList::CellType::const_iterator it0 = thisCell.begin();
                                              it0 != thisCell.end(); ++it0)
      {
        Vec3 const& xyz0 = it0->ImageCoords();
        double q0 = Charge_[it0->Idx()];
        // Exclusion list for this atom
        Iset const& excluded = Excluded_[it0->Idx()];
        // Calc interaction of atom to all other atoms in thisCell.
        for (PairList::CellType::const_iterator it1 = it0 + 1;
                                                it1 != thisCell.end(); ++it1)
        {
          Vec3 const& xyz1 = it1->ImageCoords();
          double q1 = Charge_[it1->Idx()];
          Vec3 dxyz = xyz1 - xyz0;
          double rij2 = dxyz.Magnitude2();
          // If atom excluded, calc adjustment, otherwise calc elec. energy.
          if (excluded.find( it1->Idx() ) == excluded.end())
          {
            if ( rij2 < cut2 ) {
              double rij = sqrt( rij2 );
              double qiqj = q0 * q1;
#             ifndef _OPENMP
              t_erfc_.Start();
#             endif
              //double erfc = erfc_func(ew_coeff_ * rij);
              double erfc = ERFC(ew_coeff_ * rij);
#             ifndef _OPENMP
              t_erfc_.Stop();
#             endif
              double e_elec = qiqj * erfc / rij;
              Eelec += e_elec;
              //int ta0, ta1;
              //if (atnum0 < atnum1) {
              //  ta0=atnum0; ta1=atnum1;
              //} else {
              //  ta1=atnum0; ta0=atnum1;
              //}
              //mprintf("PELEC %6i%6i%12.5f%12.5f%12.5f\n", ta0, ta1, rij, erfc, e_elec);
            }
          } else
            e_adjust += Adjust(q0, q1, sqrt(rij2));
        } // END loop over other atoms in thisCell
        // Loop over all neighbor cells
        for (unsigned int nidx = 1; nidx != cellList.size(); nidx++)
        {
          PairList::CellType const& nbrCell = PL.Cell( cellList[nidx] );
          // Translate vector for neighbor cell
          Vec3 const& tVec = PL.TransVec( transList[nidx] );
          //mprintf("\tNEIGHBOR %i (idxs %i - %i)\n", nbrCell, beg1, end1);
          // Loop over every atom in nbrCell
          for (PairList::CellType::const_iterator it1 = nbrCell.begin();
                                                  it1 != nbrCell.end(); ++it1)
          {
            Vec3 const& xyz1 = it1->ImageCoords();
            double q1 = Charge_[it1->Idx()];
            Vec3 dxyz = xyz1 + tVec - xyz0;
            double rij2 = dxyz.Magnitude2();
            //mprintf("\t\tNbrAtom %06i\n",atnum1);
            // If atom excluded, calc adjustment, otherwise calc elec. energy.
            // TODO Is there better way of checking this?
            if (excluded.find( it1->Idx() ) == excluded.end())
            {
              //mprintf("\t\t\tdist= %f\n", sqrt(rij2));
              if ( rij2 < cut2 ) {
                double rij = sqrt( rij2 );
                double qiqj = q0 * q1;
#               ifndef _OPENMP
                t_erfc_.Start();
#               endif
                //double erfc = erfc_func(ew_coeff_ * rij);
                double erfc = ERFC(ew_coeff_ * rij);
#               ifndef _OPENMP
                t_erfc_.Stop();
#               endif
                double e_elec = qiqj * erfc / rij;
                Eelec += e_elec;
                //mprintf("EELEC %4i%4i%12.5f%12.5f%12.5f%3.0f%3.0f%3.0f\n",
                //int ta0, ta1;
                //if (atnum0 < atnum1) {
                //  ta0=atnum0; ta1=atnum1;
                //} else {
                //  ta1=atnum0; ta0=atnum1;
                //}
                //mprintf("PELEC %6i%6i%12.5f%12.5f%12.5f\n", ta0, ta1, rij, erfc, e_elec);
              }
            } else
              e_adjust += Adjust(q0, q1, sqrt(rij2));
          } // END loop over neighbor cell atoms
        } // END Loop over neighbor cells
      } // Loop over thisCell atoms
    } // END if thisCell is not empty
  } // Loop over cells
# ifdef _OPENMP
  } // END pragma omp parallel
# endif
  t_direct_.Stop();
  e_adjust_out = e_adjust;
  return Eelec;
}

/** Calculate Ewald energy. Faster version that uses pair list. */
double Ewald::CalcEnergy(Frame const& frameIn, AtomMask const& maskIn)
{
  t_total_.Start();
  Matrix_3x3 ucell, recip;
  double volume = frameIn.BoxCrd().ToRecip(ucell, recip);
  double e_self = Self( volume );

  pairList_.CreatePairList(frameIn, ucell, recip, maskIn);

//  MapCoords(frameIn, ucell, recip, maskIn);
  double e_recip = Recip_Regular( recip, volume );
  double e_adjust = 0.0;
  double e_direct = Direct( pairList_, e_adjust );
  if (debug_ > 0)
    mprintf("DEBUG: Eself= %20.10f   Erecip= %20.10f   Edirect= %20.10f  Eadjust= %20.10f\n",
            e_self, e_recip, e_direct, e_adjust);
  t_total_.Stop();
  return e_self + e_recip + e_direct + e_adjust;
}

#ifdef DEBUG_EWALD
/** Calculate direct space energy. This is the slow version that doesn't
  * use a pair list; for debug purposes only.
  */
double Ewald::Direct(Matrix_3x3 const& ucell, Topology const& tIn, AtomMask const& mask)
{
  t_direct_.Start();
  double cut2 = cutoff_ * cutoff_;
  double Eelec = 0.0;
  Varray const& Image = pairList_.ImageCoords();
  Varray const& Frac = pairList_.FracCoords();
  unsigned int maxidx = Image.size();
  for (unsigned int idx1 = 0; idx1 != maxidx; idx1++)
  {
    // Set up coord for this atom
    Vec3 const& crd1 = Image[idx1];
    // Set up exclusion list for this atom
    int atom1 = mask[idx1];
    Atom::excluded_iterator excluded_atom = tIn[atom1].excludedbegin();
    for (unsigned int idx2 = idx1 + 1; idx2 != maxidx; idx2++)
    {
      int atom2 = mask[idx2];
      // If atom is excluded, just increment to next excluded atom.
      if (excluded_atom != tIn[atom1].excludedend() && atom2 == *excluded_atom) {
        ++excluded_atom;
        //mprintf("ATOM: Atom %4i to %4i excluded.\n", atom1+1, atom2+1);
      } else {
        // Only need to check nearest neighbors.
        Vec3 const& frac2 = Frac[idx2];
        for (Varray::const_iterator ixyz = Cells_.begin(); ixyz != Cells_.end(); ++ixyz)
        {
          Vec3 dxyz = ucell.TransposeMult(frac2 + *ixyz) - crd1;
          double rij2 = dxyz.Magnitude2();
          if ( rij2 < cut2 ) {
            double rij = sqrt( rij2 );
            // Coulomb
            double qiqj = Charge_[idx1] * Charge_[idx2];
            t_erfc_.Start();
            //double erfc = erfc_func(ew_coeff_ * rij);
            double erfc = ERFC(ew_coeff_ * rij);
            t_erfc_.Stop();
            double e_elec = qiqj * erfc / rij;
            Eelec += e_elec;
            //mprintf("EELEC %4i%4i%12.5f%12.5f%12.5f%3.0f%3.0f%3.0f\n",
//            mprintf("EELEC %6i%6i%12.5f%12.5f%12.5f\n", atom1, atom2, rij, erfc, e_elec);
            // TODO can we break here?
          } //else
            //mprintf("ATOM: Atom %4i to %4i outside cut, %6.2f > %6.2f %3.0f%3.0f%3.0f\n",
            //mprintf("ATOM: Atom %4i to %4i outside cut, %6.2f > %6.2f\n",
            //         atom1, atom2,sqrt(rij2),cutoff_);
        }
      }
    }
  }
  t_direct_.Stop();
  return Eelec;
}

/** Calculate Ewald energy.
  * Slow version that does not use pair list. Note that pair list is still
  * called since we require the fractional and imaged coords.
  * FIXME No Eadjust calc.
  */
double Ewald::CalcEnergy_NoPairList(Frame const& frameIn, Topology const& topIn,
                                    AtomMask const& maskIn)
{
  t_total_.Start();
  Matrix_3x3 ucell, recip;
  double volume = frameIn.BoxCrd().ToRecip(ucell, recip);
  double e_self = Self( volume );
  // Place atoms in pairlist. This calcs frac/imaged coords.
  pairList_.CreatePairList(frameIn, ucell, recip, maskIn);

  double e_recip = Recip_Regular( recip, volume );

  double e_direct = Direct( ucell, topIn, maskIn );

  //mprintf("DEBUG: Eself= %20.10f   Erecip= %20.10f   Edirect= %20.10f\n",
  //        e_self, e_recip, e_direct);
  t_total_.Stop();
  return e_self + e_recip + e_direct;
}
#endif

// Ewald::Timing()
void Ewald::Timing(double total) const {
  t_total_.WriteTiming(1,  "EwaldTotal:", total);
  t_self_.WriteTiming(2,   "Self:      ", t_total_.Total());
  t_recip_.WriteTiming(2,  "Recip:     ", t_total_.Total());
  t_trig_tables_.WriteTiming(3, "Calc trig tables:", t_recip_.Total());
  t_direct_.WriteTiming(2, "Direct :    ", t_total_.Total());
# ifndef _OPENMP
  t_erfc_.WriteTiming(3,  "ERFC:  ", t_direct_.Total());
  t_adjust_.WriteTiming(3,"Adjust:", t_direct_.Total());
# endif
  pairList_.Timing(total);
}

#include <cmath>
#include "Ewald.h"
#include "CpptrajStdio.h"
#include "Constants.h"

Ewald::Ewald() :
  sumq_(0.0),
  sumq2_(0.0),
  ew_coeff_(0.0),
  maxexp_(0.0),
  cutoff_(0.0),
  dsumTol_(0.0),
  rsumTol_(0.0),
  maxmlim_(0),
  needSumQ_(true)
{
  mlimit_[0] = 0;
  mlimit_[1] = 0;
  mlimit_[2] = 0;
  /// Save fractional translations for 1 cell in each direction (and primary cell).
  Cells_.reserve( 27 );
  for (int ix = -1; ix < 2; ix++)
    for (int iy = -1; iy < 2; iy++)
      for (int iz = -1; iz < 2; iz++)
        Cells_.push_back( Vec3(ix, iy, iz) );
}

double Ewald::INVSQRTPI_ = 1.0 / sqrt(Constants::PI);

static inline double DABS(double xIn) { if (xIn < 0.0) return -xIn; else return xIn; }
static inline int    IABS(int    xIn) { if (xIn < 0  ) return -xIn; else return xIn; }

// Original code: SANDER: erfcfun.F90
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

// Original Code: SANDER: findewaldcof
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
  mprintf("DEBUG: Ewald coefficient for cut=%g, direct sum tol=%g is %g\n",
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
  mprintf("DEBUG: MaxExp for ewcoeff=%g, direct sum tol=%g is %g\n",
          ewCoeff, rsumTol, xval);
  return xval;
}

/** Get mlimits. */
void Ewald::GetMlimits(int* mlimit, double maxexp, double eigmin, 
                       Vec3 const& reclng, Matrix_3x3 const& recip)
{
  mprintf("DEBUG: Recip lengths %12.4f%12.4f%12.4f\n", reclng[0], reclng[1], reclng[2]);

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
  mprintf("DEBUG: Number of reciprocal vectors: %i\n", nrecvecs);
}


// -----------------------------------------------------------------------------
/** Set up parameters. */
int Ewald::SetupParams(Box const& boxIn, double cutoffIn, double dsumTolIn, double rsumTolIn,
                       double ew_coeffIn, double maxexpIn,
                       const int* mlimitsIn)
{
  needSumQ_ = true;
  cutoff_ = cutoffIn;
  dsumTol_ = dsumTolIn;
  rsumTol_ = rsumTolIn;
  ew_coeff_ = ew_coeffIn;
  maxexp_ = maxexpIn;
  Matrix_3x3 ucell, recip;
  boxIn.ToRecip(ucell, recip);
  if (mlimitsIn != 0)
    std::copy(mlimitsIn, mlimitsIn+3, mlimit_);

  // Check input
  if (cutoff_ < Constants::SMALL) {
    mprinterr("Error: Direct space cutoff (%g) is too small.\n", cutoff_);
    return 1;
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

  // Set defaults if necessary
  if (DABS(ew_coeff_) < Constants::SMALL)
    ew_coeff_ = FindEwaldCoefficient( cutoff_, dsumTol_ );
  if (dsumTol_ < Constants::SMALL)
    dsumTol_ = 1E-5;
  if (rsumTol_ < Constants::SMALL)
    rsumTol_ = 5E-5;
  if (maxmlim_ > 0)
    maxexp_ = FindMaxexpFromMlim(mlimit_, recip);
  else {
    if ( maxexp_ < Constants::SMALL )
      maxexp_ = FindMaxexpFromTol(ew_coeff_, rsumTol_);
    // Calculate lengths of reciprocal vectors
    Vec3 reclng( 1.0/sqrt(recip[0]*recip[0] + recip[1]*recip[1] + recip[2]*recip[2]),
                 1.0/sqrt(recip[3]*recip[3] + recip[4]*recip[4] + recip[5]*recip[5]),
                 1.0/sqrt(recip[6]*recip[6] + recip[7]*recip[7] + recip[8]*recip[8]) );
    // eigmin typically bigger than this unless cell is badly distorted.
    double eigmin = 0.5;
    GetMlimits(mlimit_, maxexp_, eigmin, reclng, recip);
    maxmlim_ = mlimit_[0];
    maxmlim_ = std::max(maxmlim_, mlimit_[1]);
    maxmlim_ = std::max(maxmlim_, mlimit_[2]);
  }

  mprintf("DEBUG: Ewald params:\n");
  mprintf("DEBUG:   cutoff= %g   direct sum tol= %g   Ewald coeff.= %g\n",
          cutoff_, dsumTol_, ew_coeff_);
  mprintf("DEBUG:   maxexp= %g   recip. sum tol= %g\n",
          maxexp_, rsumTol_);
  mprintf("DEBUG:   mlimits= {%i,%i,%i} Max=%i\n", mlimit_[0], mlimit_[1], mlimit_[2], maxmlim_);
  return 0;
}

/** Take Cartesian coords of input atoms and map to fractional coords. */
void Ewald::MapCoords(Frame const& frmIn, Matrix_3x3 const& ucell,
                      Matrix_3x3 const& recip, AtomMask const& maskIn)
{
  Frac_.clear();
  Frac_.reserve( maskIn.Nselected() );
  Image_.clear();
  Image_.reserve( maskIn.Nselected() );

  for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom)
  {
    Vec3 fc = recip * Vec3(frmIn.XYZ(*atom));
    // Wrap back into primary cell
    //Frac_.push_back( Vec3(fc[0]-floor(fc[0]), fc[1]-floor(fc[1]), fc[2]-floor(fc[2])) );
    Frac_.push_back( Vec3(fc[0]-(int)fc[0], fc[1]-(int)fc[1], fc[2]-(int)fc[2]) );
    Image_.push_back( ucell.TransposeMult( Frac_.back() ) );
  }
  mprintf("DEBUG: Mapped coords for %zu atoms.\n", Frac_.size());
}

/** Calculate sum of charges and squared charges. */
void Ewald::CalcSumQ(Topology const& topIn, AtomMask const& maskIn) {
  sumq_ = 0.0;
  sumq2_ = 0.0;
  Charge_.clear();
  for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom) {
    double qi = topIn[*atom].Charge() * Constants::ELECTOAMBER;
    Charge_.push_back(qi);
    sumq_ += qi;
    sumq2_ += (qi * qi);
  }
  mprintf("DEBUG: sumq= %20.10f   sumq2= %20.10f\n", sumq_, sumq2_);
  needSumQ_ = false;
}

/** Self energy. */
double Ewald::Self(double volume) {
  double d0 = -ew_coeff_ * INVSQRTPI_;
  double ene = sumq2_ * d0;
  mprintf("DEBUG: d0= %20.10f   ene= %20.10f\n", d0, ene);
  double factor = Constants::PI / (ew_coeff_ * ew_coeff_ * volume);
  double ee_plasma = -0.5 * factor * sumq_ * sumq_;
  ene += ee_plasma;
  return ene;
}

/** Recip energy. */
double Ewald::Recip_Regular(Matrix_3x3 const& recip, double volume) {
  double fac = (Constants::PI*Constants::PI) / (ew_coeff_ * ew_coeff_);
  double maxexp2 = maxexp_ * maxexp_;
  double ene = 0.0;
  // Number of M values is the max + 1.
  int mmax = maxmlim_ + 1;
  // Build exponential factors for use in structure factors.
  // These arrays are laid out in 1D; value for each atom at each m, i.e.
  // A0M0 A1M0 A2M0 ... ANM0 A0M1 ... ANMX
  Darray cosf1, cosf2, cosf3, sinf1, sinf2, sinf3;
  cosf1.reserve( Frac_.size()*mmax );
  cosf2.reserve( Frac_.size()*mmax );
  cosf3.reserve( Frac_.size()*mmax );
  sinf1.reserve( Frac_.size()*mmax );
  sinf2.reserve( Frac_.size()*mmax );
  sinf3.reserve( Frac_.size()*mmax );
  // M0
  for (unsigned int i = 0; i != Frac_.size(); i++) {
    cosf1.push_back( 1.0 );
    cosf2.push_back( 1.0 );
    cosf3.push_back( 1.0 );
    sinf1.push_back( 0.0 );
    sinf2.push_back( 0.0 );
    sinf3.push_back( 0.0 );
  }
  // M1
  for (unsigned int i = 0; i != Frac_.size(); i++) {
    //mprintf("FRAC: %6i%20.10f%20.10f%20.10f\n", i+1, Frac_[i][0], Frac_[i][1], Frac_[i][2]);
    cosf1.push_back( cos(Constants::TWOPI * Frac_[i][0]) );
    cosf2.push_back( cos(Constants::TWOPI * Frac_[i][1]) );
    cosf3.push_back( cos(Constants::TWOPI * Frac_[i][2]) );
    sinf1.push_back( sin(Constants::TWOPI * Frac_[i][0]) );
    sinf2.push_back( sin(Constants::TWOPI * Frac_[i][1]) );
    sinf3.push_back( sin(Constants::TWOPI * Frac_[i][2]) );
  }
  // M2-MX
  // Get the higher factors by recursion using trig addition rules.
  // Negative values of M by complex conjugation, or even cosf, odd sinf.
  // idx will always point to M-1 values
  unsigned int idx = Frac_.size();
  for (int m = 2; m < mmax; m++) {
    // Set m1idx to beginning of M1 values.
    unsigned int m1idx = Frac_.size();
    for (unsigned int i = 0; i != Frac_.size(); i++, idx++, m1idx++) {
      cosf1.push_back( cosf1[idx]*cosf1[m1idx] - sinf1[idx]*sinf1[m1idx] );
      cosf2.push_back( cosf2[idx]*cosf2[m1idx] - sinf2[idx]*sinf2[m1idx] );
      cosf3.push_back( cosf3[idx]*cosf3[m1idx] - sinf3[idx]*sinf3[m1idx] );
      sinf1.push_back( sinf1[idx]*cosf1[m1idx] + cosf1[idx]*sinf1[m1idx] );
      sinf2.push_back( sinf2[idx]*cosf2[m1idx] + cosf2[idx]*sinf2[m1idx] );
      sinf3.push_back( sinf3[idx]*cosf3[m1idx] + cosf3[idx]*sinf3[m1idx] );
    }
  }
  // DEBUG
/*  unsigned int midx = 0;
  for (int m = 0; m != mmax; m++) {
    for (unsigned int i = 0; i != Frac_.size(); i++, midx++)
      mprintf("TRIG: %6i%6u%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n", m,i+1,
               cosf1[midx], cosf2[midx], cosf3[midx],
               sinf1[midx], sinf2[midx], sinf3[midx]);
  }*/

  double mult = 1.0;
//  int count = -1;
  Darray c12(Frac_.size(), 0.0);
  Darray s12(Frac_.size(), 0.0);
  Darray c3(Frac_.size(), 0.0);
  Darray s3(Frac_.size(), 0.0);
  for (int m1 = 0; m1 <= mlimit_[0]; m1++)
  {
    for (int m2 = -mlimit_[1]; m2 <= mlimit_[1]; m2++)
    {
//      count++;
//      if ( iproc == (count % nproc) )
//      {
        int m1idx = Frac_.size() * m1;
        int m2idx = Frac_.size() * IABS(m2);
        if (m2 < 0) {
          for (unsigned int i = 0; i != Frac_.size(); i++, m1idx++, m2idx++) {
            c12[i] = cosf1[m1idx]*cosf2[m2idx] + sinf1[m1idx]*sinf2[m2idx];
            s12[i] = sinf1[m1idx]*cosf2[m2idx] - cosf1[m1idx]*sinf2[m2idx];
          }
        } else {
          for (unsigned int i = 0; i != Frac_.size(); i++, m1idx++, m2idx++) {
            c12[i] = cosf1[m1idx]*cosf2[m2idx] - sinf1[m1idx]*sinf2[m2idx];
            s12[i] = sinf1[m1idx]*cosf2[m2idx] + cosf1[m1idx]*sinf2[m2idx];
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
//          double vterm = 0.0;
          if ( m1*m1 + m2*m2 + m3*m3 > 0 ) {
            eterm = exp(-fac*msq) / denom;
//            vterm = 2.0 * (fac*msq + 1.0) / msq;
          }
          // mult takes care to double count for symmetry. Can take care of
          // with eterm.
          eterm *= mult;
          if (msq < maxexp2) {
            int m3idx = Frac_.size() * IABS(m3);
            // Get the product of complex exponentials.
            if (m3 < 0) {
              for (unsigned int i = 0; i != Frac_.size(); i++, m3idx++) {
                c3[i] = c12[i]*cosf3[m3idx] + s12[i]*sinf3[m3idx];
                s3[i] = s12[i]*cosf3[m3idx] - c12[i]*sinf3[m3idx];
              }
            } else {
              for (unsigned int i = 0; i != Frac_.size(); i++, m3idx++) {
                c3[i] = c12[i]*cosf3[m3idx] - s12[i]*sinf3[m3idx];
                s3[i] = s12[i]*cosf3[m3idx] + c12[i]*sinf3[m3idx];
              }
            }
            // Get the structure factor
            double cstruct = 0.0;
            double sstruct = 0.0;
            for (unsigned int i = 0; i != Frac_.size(); i++) {
              cstruct += Charge_[i] * c3[i];
              sstruct += Charge_[i] * s3[i];
            }
            double struc2 = cstruct*cstruct + sstruct*sstruct;
            ene += eterm * struc2;
            //mprintf("LOOP: %3i%3i%3i ENE= %20.10f\n", m1, m2, m3, ene);
          } // END IF msq < maxexp2
        } // END loop over m3
//      } // END IF iproc == (count % nproc)
    } // END loop over m2
    mult = 2.0;
  } // END loop over m1

  return ene * 0.5;
}

/** Direct space energy. */
double Ewald::Direct(Matrix_3x3 const& ucell, Topology const& tIn, AtomMask const& mask)
{
//  time_NB_.Start();
  double cut2 = cutoff_ * cutoff_;
  double Eelec = 0.0;
  unsigned int maxidx = Image_.size();
  for (unsigned int idx1 = 0; idx1 != maxidx; idx1++)
  {
    // Set up coord for this atom
    Vec3 const& crd1 = Image_[idx1];
    // Set up exclusion list for this atom
    int atom1 = mask[idx1];
    Atom::excluded_iterator excluded_atom = tIn[atom1].excludedbegin();
    for (unsigned int idx2 = idx1 + 1; idx2 != maxidx; idx2++)
    {
      int atom2 = mask[idx2];
      // If atom is excluded, just increment to next excluded atom.
      mprintf("ATOM: Atom %4i to %4i", atom1+1, atom2+1);
      if (excluded_atom != tIn[atom1].excludedend() && atom2 == *excluded_atom) {
        ++excluded_atom;
        mprintf(" excluded.\n");
      } else {
        // Only need to check nearest neighbors.
        Vec3 const& frac2 = Frac_[idx2];
        for (Varray::const_iterator ixyz = Cells_.begin(); ixyz != Cells_.end(); ++ixyz)
        {
          Vec3 dxyz = ucell.TransposeMult(frac2 + *ixyz) - crd1;
          double rij2 = dxyz.Magnitude2();
          if ( rij2 < cut2 ) {
            mprintf("\n");
            double rij = sqrt( rij2 );
            // Coulomb
            double qiqj = Charge_[idx1] * Charge_[idx2];
            double erfc = erfc_func(ew_coeff_ * rij);
            double e_elec = qiqj * erfc / rij;
            Eelec += e_elec;
            mprintf("EELEC %4i%4i%12.5f%12.5f%12.5f\n", atom1+1, atom2+1, rij, erfc, e_elec);
          } else
            mprintf(" outside cut, %g > %g.\n", sqrt(rij2), cutoff_);
        }
      }
    }
  }
//  time_NB_.Stop();
  return Eelec;
}

/** Calculate Ewald energy. */
double Ewald::CalcEnergy(Frame const& frameIn, Topology const& topIn, AtomMask const& maskIn)
{
  Matrix_3x3 ucell, recip;
  double volume = frameIn.BoxCrd().ToRecip(ucell, recip);
  double e_self = Self( volume );

  MapCoords(frameIn, ucell, recip, maskIn);
  double e_recip = Recip_Regular( recip, volume );

  double e_direct = Direct( ucell, topIn, maskIn );

  mprintf("DEBUG: Eself= %20.10f   Erecip= %20.10f   Edirect= %20.10f\n",
          e_self, e_recip, e_direct);

  return e_self + e_recip + e_direct;
}

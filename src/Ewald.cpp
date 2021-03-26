#include <cmath> //sqrt
#include "Ewald.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "StringRoutines.h" // ByteString
#include "Spline.h"
#include "Topology.h"
#include "CharMask.h"
#include "ParameterTypes.h"
#ifdef DEBUG_PAIRLIST
#incl ude "PDBfile.h"
#endif

/// CONSTRUCTOR
Ewald::Ewald() :
  ew_coeff_(0.0),
  lw_coeff_(0.0),
  switch_width_(0.0),
  cutoff_(0.0),
  cut2_(0.0),
  cut2_0_(0.0),
  dsumTol_(0.0),
  debug_(0),
  sumq_(0.0),
  sumq2_(0.0),
  Vdw_Recip_term_(0)
{
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

// Ewald::ERFC()
double Ewald::ERFC(double xIn) const {
  return table_.Yval( xIn);
}

/** Non-inlined version of ERFC */
double Ewald::ErfcFxn(double xIn) const {
  return table_.Yval( xIn );
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

/** Convert charges to Amber units. Calculate sum of charges and squared charges. */
void Ewald::CalculateCharges(Topology const& topIn, AtomMask const& maskIn) {
  sumq_ = 0.0;
  sumq2_ = 0.0;
  Charge_.clear();
  TypeIndices_.clear();
  for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom) {
    double qi = topIn[*atom].Charge() * Constants::ELECTOAMBER;
    Charge_.push_back(qi);
    sumq_ += qi;
    sumq2_ += (qi * qi);
    // Store atom type indices for selected atoms.
    TypeIndices_.push_back( topIn[*atom].TypeIndex() );
  }
  //mprintf("DEBUG: sumq= %20.10f   sumq2= %20.10f\n", sumq_, sumq2_);
  Setup_VDW_Correction( topIn, maskIn );
}

void Ewald::CalculateC6params(Topology const& topIn, AtomMask const& maskIn) {
  Cparam_.clear();
  if (lw_coeff_ > 0.0) {
    for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom)
    {
      double rmin = topIn.GetVDWradius( *atom );
      double eps  = topIn.GetVDWdepth( *atom );
      Cparam_.push_back( 8.0 * (rmin*rmin*rmin) * sqrt(2 * eps) );
      if (debug_ > 0)
        mprintf("DEBUG: C6 param atom %8i = %16.8f\n", *atom+1, Cparam_.back());
    }
  } else
    Cparam_.assign(maskIn.Nselected(), 0.0);
}

/** Set up exclusion lists for selected atoms. */
void Ewald::SetupExclusionList(Topology const& topIn, AtomMask const& maskIn)
{
  // Use distance of 4 (up to dihedrals)
  if (Excluded_.SetupExcluded(topIn.Atoms(), maskIn, 4,
                              ExclusionArray::EXCLUDE_SELF,
                              ExclusionArray::FULL))
  {
    mprinterr("Error: Ewald: Could not set up exclusion list.\n");
    return;
  }
}

/** Check some common input. */
int Ewald::CheckInput(Box const& boxIn, int debugIn, double cutoffIn, double dsumTolIn,
                      double ew_coeffIn, double lw_coeffIn, double switch_widthIn,
                      double erfcTableDxIn, double skinnbIn)
{
  debug_ = debugIn;
  cutoff_ = cutoffIn;
  dsumTol_ = dsumTolIn;
  ew_coeff_ = ew_coeffIn;
  lw_coeff_ = lw_coeffIn;
  switch_width_ = switch_widthIn;
  double erfcTableDx = erfcTableDxIn;
  // Check input
  if (cutoff_ < Constants::SMALL) {
    mprinterr("Error: Direct space cutoff (%g) is too small.\n", cutoff_);
    return 1;
  }
  char dir[3] = {'X', 'Y', 'Z'};
  // NOTE: First 3 box parameters are X Y Z
  for (int i = 0; i < 3; i++) {
    if (cutoff_ > boxIn.Param((Box::ParamType)i)/2.0) {
      mprinterr("Error: Cutoff must be less than half the box length (%g > %g, %c)\n",
                cutoff_, boxIn.Param((Box::ParamType)i)/2.0, dir[i]);
      return 1;
    }
  }
  if (skinnbIn < 0.0) {
    mprinterr("Error: skinnb is less than 0.0\n");
    return 1;
  }
  if (switch_width_ < 0.0) switch_width_ = 0.0;
  if (switch_width_ > cutoff_) {
    mprinterr("Error: Switch width must be less than the cutoff.\n");
    return 1;
  }

  // Set defaults if necessary
  if (dsumTol_ < Constants::SMALL)
    dsumTol_ = 1E-5;
  if (DABS(ew_coeff_) < Constants::SMALL)
    ew_coeff_ = FindEwaldCoefficient( cutoff_, dsumTol_ );
  if (erfcTableDx <= 0.0) erfcTableDx = 1.0 / 5000;
  // TODO make this optional
  if (table_.FillTable( erfc_func, erfcTableDx, 0.0, cutoff_*ew_coeff_*1.5 )) {
    mprinterr("Error: Could not set up spline table for ERFC\n");
    return 1;
  }
  table_.PrintMemUsage("\t");
  table_.PrintTableInfo("\t");
  // TODO do for C6 as well
  // TODO for C6 correction term
  if (lw_coeff_ < 0.0)
    lw_coeff_ = 0.0;
  else if (DABS(lw_coeff_) < Constants::SMALL)
    lw_coeff_ = ew_coeff_;

  // Calculate some common factors.
  cut2_ = cutoff_ * cutoff_;
  double cut0 = cutoff_ - switch_width_;
  cut2_0_ = cut0 * cut0;

  return 0;
}

/** Initialize and set up pairlist. */
int Ewald::Setup_Pairlist(Box const& boxIn, double skinnbIn) {
  if (pairList_.InitPairList(cutoff_, skinnbIn, debug_)) return 1;
  if (pairList_.SetupPairList( boxIn )) return 1;
# ifdef DEBUG_PAIRLIST
  // Write grid PDB
  PDBfile gridpdb;
  gridpdb.OpenWrite("gridpoints.pdb");
  for (int iz = 0; iz != pairList_.NZ(); iz++)
    for (int iy = 0; iy != pairList_.NY(); iy++)
      for (int ix = 0; ix != pairList_.NX(); ix++) {
        double fx = (double)ix / (double)pairList_.NX();
        double fy = (double)iy / (double)pairList_.NY();
        double fz = (double)iz / (double)pairList_.NZ();
        Vec3 cart = boxIn.UnitCell().TransposeMult( Vec3(fx,fy,fz) );
        gridpdb.WriteHET(1, cart[0], cart[1], cart[2]);
      }
  gridpdb.CloseFile();
# endif
  return 0;
}

/** Electrostatic self energy. This is the cancelling Gaussian plus the "neutralizing plasma". */
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

/** Lennard-Jones self energy. */
double Ewald::Self6() {
  t_self_.Start(); // TODO precalc
  double ew2 = lw_coeff_ * lw_coeff_;
  double ew6 = ew2 * ew2 * ew2;
  double c6sum = 0.0;
  for (Darray::const_iterator it = Cparam_.begin(); it != Cparam_.end(); ++it)
    c6sum += ew6 * (*it * *it);
  t_self_.Stop();
  return c6sum / 12.0;
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

/** Ewald adjustment, for inheriting classes. */
#ifdef _OPENMP
double Ewald::AdjustFxn(double q0, double q1, double rij) const {
  return Adjust(q0, q1, rij);
}
#else
double Ewald::AdjustFxn(double q0, double q1, double rij) {
  return Adjust(q0, q1, rij);
}
#endif

/** Switching function for Lennard-Jones. */
static inline double switch_fn(double rij2, double cut2_0, double cut2_1)
{
  if (rij2 <= cut2_0)
    return 1.0;
  else if (rij2 > cut2_1)
    return 0.0;
  else {
    double xoff_m_x = cut2_1 - rij2;
    double fac = 1.0 / (cut2_1 - cut2_0);
    return (xoff_m_x*xoff_m_x) * (cut2_1 + 2.0*rij2 - 3.0*cut2_0) * (fac*fac*fac);
  }
}

/** Switching function for Lennard-Jones, for inheriting classes. */
double Ewald::SwitchFxn(double rij2, double cut2_0, double cut2_1) {
  return switch_fn(rij2, cut2_0, cut2_1);
}

/** Nonbond direct-space calculation for Coulomb electrostatics and Lennard-Jones,
  * intended for use with long-range LJ correction.
  */
double Ewald::Direct_VDW_LongRangeCorrection(PairList const& PL, double& evdw_out)
{
  t_direct_.Start();
  double Eelec = 0.0;
  double e_adjust = 0.0;
  double Evdw = 0.0;
  int cidx;
# ifdef _OPENMP
# pragma omp parallel private(cidx) reduction(+: Eelec, Evdw, e_adjust)
  {
# pragma omp for
# endif
# include "PairListLoop.h"
# ifdef _OPENMP
  } // END pragma omp parallel
# endif
  t_direct_.Stop();
# ifdef DEBUG_PAIRLIST
  mprintf("DEBUG: Elec                             = %16.8f\n", Eelec);
  mprintf("DEBUG: Eadjust                          = %16.8f\n", e_adjust);
  mprintf("DEBUG: LJ vdw                           = %16.8f\n", Evdw);
# endif
  evdw_out = Evdw;
  return Eelec + e_adjust;
}

/** Nonbond direct-space calculation for Coulomb electrostatics and Lennard-Jones
  * calculated via PME.
  */
double Ewald::Direct_VDW_LJPME(PairList const& PL, double& evdw_out)
{
  t_direct_.Start();
  double Eelec = 0.0;
  double e_adjust = 0.0;
  double Evdw = 0.0;
  double Eljpme_correction = 0.0;
  double Eljpme_correction_excl = 0.0;
  int cidx;
# define CPPTRAJ_EKERNEL_LJPME
# ifdef _OPENMP
# pragma omp parallel private(cidx) reduction(+: Eelec, Evdw, e_adjust, Eljpme_correction,Eljpme_correction_excl)
  {
# pragma omp for
# endif
# include "PairListLoop.h"
# ifdef _OPENMP
  } // END pragma omp parallel
# endif
# undef CPPTRAJ_EKERNEL_LJPME
  t_direct_.Stop();
# ifdef DEBUG_PAIRLIST
  mprintf("DEBUG: Elec                             = %16.8f\n", Eelec);
  mprintf("DEBUG: Eadjust                          = %16.8f\n", e_adjust);
  mprintf("DEBUG: LJ vdw                           = %16.8f\n", Evdw);
  mprintf("DEBUG: LJ vdw PME correction            = %16.8f\n", Eljpme_correction);
  mprintf("DEBUG: LJ vdw PME correction (excluded) = %16.8f\n", Eljpme_correction_excl);
# endif
  evdw_out = Evdw + Eljpme_correction + Eljpme_correction_excl;
  return Eelec + e_adjust;
}


//  Ewald::Direct()
/** Calculate direct space energy. This is the faster version that uses
  * a pair list. Also calculate the energy adjustment for excluded
  * atoms.
  * \param PL The pairlist used to calculate energy.
  * \param e_adjust_out The electrostatic adjust energy for excluded atoms.
  * \param evdw_out The direct space van der Waals term (corrected for exclusion if LJ PME).
  * \return The electrostatics term plus exclusion adjustment.
  */
double Ewald::Direct(PairList const& PL, double& evdw_out)
{
  if (lw_coeff_ > 0.0)
    return Direct_VDW_LJPME(PL, evdw_out);
  else
    return Direct_VDW_LongRangeCorrection(PL, evdw_out);
}

/** Determine VDW long range correction prefactor. */
void Ewald::Setup_VDW_Correction(Topology const& topIn, AtomMask const& maskIn) {
  Vdw_Recip_term_ = 0.0;
  NB_ = static_cast<NonbondParmType const*>( &(topIn.Nonbond()) );
  if (!NB_->HasNonbond()) {
    mprintf("Warning: '%s' has no nonbonded parameters. Cannot calculate VDW correction.\n",
            topIn.c_str());
    return;
  }
  // Count the number of each unique nonbonded type.
  N_vdw_type_.assign( NB_->Ntypes(), 0 );
  vdw_type_.clear();
  for (AtomMask::const_iterator atm = maskIn.begin(); atm != maskIn.end(); ++atm)
  {
    N_vdw_type_[ topIn[*atm].TypeIndex() ]++;
    vdw_type_.push_back( topIn[*atm].TypeIndex() );
  }
  if (debug_ > 0) {
    mprintf("DEBUG: %zu VDW types.\n", N_vdw_type_.size());
    for (Iarray::const_iterator it = N_vdw_type_.begin(); it != N_vdw_type_.end(); ++it)
      mprintf("\tType %li = %i\n", it-N_vdw_type_.begin(), *it);
  }
  // Determine correction term from types and LJ B parameters
  for (unsigned int itype = 0; itype != N_vdw_type_.size(); itype++)
  {
    double atype_vdw_term = 0.0; // term for each nonbond atom type
    unsigned int offset = N_vdw_type_.size() * itype;
    for (unsigned int jtype = 0; jtype != N_vdw_type_.size(); jtype++)
    {
      unsigned int idx = offset + jtype;
      int nbidx = NB_->NBindex()[ idx ];
      if (nbidx > -1) {
        atype_vdw_term += N_vdw_type_[itype] * N_vdw_type_[jtype] * NB_->NBarray()[ nbidx ].B();

        Vdw_Recip_term_ += N_vdw_type_[itype] * N_vdw_type_[jtype] * NB_->NBarray()[ nbidx ].B();
      }
    }
    atype_vdw_recip_terms_.push_back(atype_vdw_term);  // the nonbond interaction for each atom type
  }
}

/** Calculate full VDW long range correction from volume. */
double Ewald::Vdw_Correction(double volume) {
  double prefac = Constants::TWOPI / (3.0*volume*cutoff_*cutoff_*cutoff_);
  double e_vdwr = -prefac * Vdw_Recip_term_;
  if (debug_ > 0) mprintf("DEBUG: Vdw correction %20.10f\n", e_vdwr);
  return e_vdwr;
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
  double volume = frameIn.BoxCrd().CellVolume();
  double e_self = Self( volume );
  // Place atoms in pairlist. This calcs frac/imaged coords.
  pairList_.CreatePairList(frameIn, frameIn.BoxCrd().UnitCell(), frameIn.BoxCrd().FracCell(), maskIn);

  double e_recip = Recip_Regular( recip, volume );

  double e_direct = Direct( ucell, topIn, maskIn );

  //mprintf("DEBUG: Eself= %20.10f   Erecip= %20.10f   Edirect= %20.10f\n",
  //        e_self, e_recip, e_direct);
  t_total_.Stop();
  return e_self + e_recip + e_direct;
}
#endif /* DEBUG_EWALD */

// Ewald::Timing()
void Ewald::Timing(double total) const {
  t_total_.WriteTiming(1,  "  EwaldTotal:", total);
  t_self_.WriteTiming(2,   "Self:      ", t_total_.Total());
  t_recip_.WriteTiming(2,  "Recip:     ", t_total_.Total());
  if (t_trig_tables_.Total() > 0.0)
    t_trig_tables_.WriteTiming(3, "Calc trig tables:", t_recip_.Total());
  t_direct_.WriteTiming(2, "Direct:    ", t_total_.Total());
# ifndef _OPENMP
  t_erfc_.WriteTiming(3,  "ERFC:  ", t_direct_.Total());
  t_adjust_.WriteTiming(3,"Adjust:", t_direct_.Total());
# endif
  pairList_.Timing(total);
}

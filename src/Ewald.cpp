#include <cmath> //sqrt
#include "Ewald.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "StringRoutines.h" // ByteString
#include "Spline.h"
#ifdef DEBUG_PAIRLIST
#incl ude "PDBfile.h"
#endif

/// CONSTRUCTOR
Ewald::Ewald() :
  sumq_(0.0),
  sumq2_(0.0),
  ew_coeff_(0.0),
  cutoff_(0.0),
  dsumTol_(0.0),
  erfcTableDx_(0.0),
  one_over_Dx_(0.0),
  debug_(0)
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

/** Set up exclusion lists for selected atoms. */
void Ewald::SetupExcluded(Topology const& topIn, AtomMask const& maskIn)
{
  Excluded_.clear();
  Excluded_.resize( maskIn.Nselected() );
  // Create a character mask so we can see if atoms in excluded lists are
  // also selected.
  CharMask Cmask(maskIn.ConvertToCharMask(), maskIn.Nselected());
  // Create a map of atom number to maskIn index.
  int selectedIdx = 0;
  Iarray atToIdx( Cmask.Natom(), -1 );
  for (int cidx = 0; cidx != Cmask.Natom(); cidx++)
    if (Cmask.AtomInCharMask(cidx))
      atToIdx[cidx] = selectedIdx++;
  // Loop over selected atoms
  for (int idx = 0; idx != maskIn.Nselected(); idx++)
  {
    // Always exclude self
    Excluded_[idx].insert( idx );
    int at = maskIn[idx];
    for (Atom::excluded_iterator excluded_atom = topIn[at].excludedbegin();
                                 excluded_atom != topIn[at].excludedend();
                               ++excluded_atom)
    {
      if (Cmask.AtomInCharMask(*excluded_atom))
      {
        // Find excluded atoms index in maskIn
        int excluded_idx = atToIdx[*excluded_atom];
        Excluded_[idx         ].insert( excluded_idx );
        Excluded_[excluded_idx].insert( idx          );
      }
    }
  }
  unsigned int ex_size = 0;
  for (Iarray2D::const_iterator it = Excluded_.begin(); it != Excluded_.end(); ++it)
    ex_size += it->size();
  mprintf("\tMemory used by full exclusion list: %s\n",
          ByteString(ex_size * sizeof(int), BYTE_DECIMAL).c_str());
}

/** Check some common input. */
int Ewald::CheckInput(Box const& boxIn, int debugIn, double cutoffIn, double dsumTolIn,
                      double ew_coeffIn, double erfcTableDxIn, double skinnbIn)
{
  debug_ = debugIn;
  cutoff_ = cutoffIn;
  dsumTol_ = dsumTolIn;
  ew_coeff_ = ew_coeffIn;
  erfcTableDx_ = erfcTableDxIn;
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
  if (skinnbIn < 0.0) {
    mprinterr("Error: skinnb is less than 0.0\n");
    return 1;
  }

  // Set defaults if necessary
  if (dsumTol_ < Constants::SMALL)
    dsumTol_ = 1E-5;
  if (DABS(ew_coeff_) < Constants::SMALL)
    ew_coeff_ = FindEwaldCoefficient( cutoff_, dsumTol_ );
  if (erfcTableDx_ <= 0.0) erfcTableDx_ = 1.0 / 5000;
  // TODO make this optional
  FillErfcTable( cutoff_, ew_coeff_ );

  return 0;
}

/** Initialize and set up pairlist. */
int Ewald::Setup_Pairlist(Box const& boxIn, Vec3 const& recipLengths, double skinnbIn) {
  if (pairList_.InitPairList(cutoff_, skinnbIn, debug_)) return 1;
  if (pairList_.SetupPairList( boxIn.Type(), recipLengths )) return 1;
# ifdef DEBUG_PAIRLIST
  // Write grid PDB
  PDBfile gridpdb;
  Matrix_3x3 ucell, recip;
  boxIn.ToRecip(ucell, recip);
  gridpdb.OpenWrite("gridpoints.pdb");
  for (int iz = 0; iz != pairList_.NZ(); iz++)
    for (int iy = 0; iy != pairList_.NY(); iy++)
      for (int ix = 0; ix != pairList_.NX(); ix++) {
        double fx = (double)ix / (double)pairList_.NX();
        double fy = (double)iy / (double)pairList_.NY();
        double fz = (double)iz / (double)pairList_.NZ();
        Vec3 cart = ucell.TransposeMult( Vec3(fx,fy,fz) );
        gridpdb.WriteHET(1, cart[0], cart[1], cart[2]);
      }
  gridpdb.CloseFile();
# endif
  return 0;
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
double Ewald::Direct(PairList const& PL, double& e_adjust_out, double& evdw_out)
{
  t_direct_.Start();
  double cut2 = cutoff_ * cutoff_;
  double Eelec = 0.0;
  double e_adjust = 0.0;
  double Evdw = 0.0;
  int cidx;
# ifdef _OPENMP
# pragma omp parallel private(cidx) reduction(+: Eelec, Evdw, e_adjust)
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
#       ifdef DEBUG_PAIRLIST
        mprintf("DBG: Cell %6i (%6i atoms):\n", cidx+1, thisCell.NatomsInGrid());
#       endif
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
#         ifdef DEBUG_PAIRLIST
          mprintf("\tAtom %6i to atom %6i (%f)\n", it0->Idx()+1, it1->Idx()+1, sqrt(rij2));
#         endif
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
              //if (it0->Idx() < it1->Idx()) {
              //  ta0=it0->Idx(); ta1=it1->Idx();
              //} else {
              //  ta1=it0->Idx(); ta0=it1->Idx();
              //}
              //mprintf("PELEC %6i%6i%12.5f%12.5f%12.5f\n", ta0, ta1, rij, erfc, e_elec);
              int nbindex = NB_->GetLJindex(TypeIndices_[it0->Idx()],
                                            TypeIndices_[it1->Idx()]);
              if (nbindex > -1) {
                NonbondType const& LJ = NB_->NBarray()[ nbindex ];
                double r2    = 1.0 / rij2;
                double r6    = r2 * r2 * r2;
                double r12   = r6 * r6;
                double f12   = LJ.A() * r12;  // A/r^12
                double f6    = LJ.B() * r6;   // B/r^6
                double e_vdw = f12 - f6;      // (A/r^12)-(B/r^6)
                Evdw += e_vdw;
                //mprintf("PVDW %8i%8i%20.6f%20.6f\n", ta0+1, ta1+1, e_vdw, r2);
              }
            }
          } else
            e_adjust += Adjust(q0, q1, sqrt(rij2));
        } // END loop over other atoms in thisCell
        // Loop over all neighbor cells
        for (unsigned int nidx = 1; nidx != cellList.size(); nidx++)
        {
          PairList::CellType const& nbrCell = PL.Cell( cellList[nidx] );
#         ifdef DEBUG_PAIRLIST
          if (nbrCell.NatomsInGrid()>0) mprintf("\tto neighbor cell %6i\n", cellList[nidx]+1);
#         endif
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
#           ifdef DEBUG_PAIRLIST
            mprintf("\t\tAtom %6i to atom %6i (%f)\n", it0->Idx()+1, it1->Idx()+1, sqrt(rij2));
#           endif
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
                //if (it0->Idx() < it1->Idx()) {
                //  ta0=it0->Idx(); ta1=it1->Idx();
                //} else {
                //  ta1=it0->Idx(); ta0=it1->Idx();
                //}
                //mprintf("PELEC %6i%6i%12.5f%12.5f%12.5f\n", ta0, ta1, rij, erfc, e_elec);
                int nbindex = NB_->GetLJindex(TypeIndices_[it0->Idx()],
                                              TypeIndices_[it1->Idx()]);
                if (nbindex > -1) {
                  NonbondType const& LJ = NB_->NBarray()[ nbindex ];
                  double r2    = 1.0 / rij2;
                  double r6    = r2 * r2 * r2;
                  double r12   = r6 * r6;
                  double f12   = LJ.A() * r12;  // A/r^12
                  double f6    = LJ.B() * r6;   // B/r^6
                  double e_vdw = f12 - f6;      // (A/r^12)-(B/r^6)
                  Evdw += e_vdw;
                  //mprintf("PVDW %8i%8i%20.6f%20.6f\n", ta0+1, ta1+1, e_vdw, r2);
                }
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
  evdw_out = Evdw;
  return Eelec;
}

/** Determine VDW long range correction prefactor. */
void Ewald::Setup_VDW_Correction(Topology const& topIn, AtomMask const& maskIn) {
  Vdw_Recip_term_ = 0.0;
  NB_ = static_cast<NonbondParmType const*>( &(topIn.Nonbond()) );
  if (!NB_->HasNonbond()) {
    mprintf("Warning: '%s' has no nonbonded parameters. Cannot calculate VDW correction.\n");
    return;
  }
  // Count the number of each unique nonbonded type.
  Iarray N_vdw_type( NB_->Ntypes(), 0 );
  for (AtomMask::const_iterator atm = maskIn.begin(); atm != maskIn.end(); ++atm)
    N_vdw_type[ topIn[*atm].TypeIndex() ]++;
  if (debug_ > 0) {
    mprintf("DEBUG: %zu VDW types.\n", N_vdw_type.size());
    for (Iarray::const_iterator it = N_vdw_type.begin(); it != N_vdw_type.end(); ++it)
      mprintf("\tType %u = %i\n", it-N_vdw_type.begin(), *it);
  }
  // Determine correction term from types and LJ B parameters
  for (unsigned int itype = 0; itype != N_vdw_type.size(); itype++)
  {
    unsigned int offset = N_vdw_type.size() * itype;
    for (unsigned int jtype = 0; jtype != N_vdw_type.size(); jtype++)
    {
      unsigned int idx = offset + jtype;
      int nbidx = NB_->NBindex()[ idx ];
      if (nbidx > -1)
        Vdw_Recip_term_ += N_vdw_type[itype] * N_vdw_type[jtype] * NB_->NBarray()[ nbidx ].B();
    }
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

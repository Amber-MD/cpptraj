#include <cmath>
#include "Energy.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "TorsionRoutines.h"
#include "Constants.h"
#include "Topology.h"
#include "CharMask.h"
#include "ExclusionArray.h"

const double Energy_Amber::QFAC = Constants::ELECTOAMBER * Constants::ELECTOAMBER;

// CONSTRUCTOR
Energy_Amber::Energy_Amber() : debug_(0) {}

/** Bond energy */
double Energy_Amber::E_bond(Frame const& fIn, Topology const& tIn, CharMask const& mask)
{
  // Heavy atom bonds
  double Ebond = CalcBondEnergy(fIn, tIn.Bonds(), tIn.BondParm(), mask);
  Ebond += CalcBondEnergy(fIn, tIn.BondsH(), tIn.BondParm(), mask);
  return Ebond;
}

// Energy_Amber::CalcBondEnergy()
double Energy_Amber::CalcBondEnergy(Frame const& fIn, BondArray const& Bonds,
                                    BondParmArray const& BPA, CharMask const& mask)
{
  double Ebond = 0.0;
  for (BondArray::const_iterator b = Bonds.begin(); b != Bonds.end(); ++b)
  {
    if (mask.AtomInCharMask( b->A1() ) && mask.AtomInCharMask( b->A2() ))
    {
      int bpidx = b->Idx();
      if (bpidx < 0) {
        if (debug_ > 0)
          mprintf("Warning: Bond %i -- %i has no parameters.\n", b->A1()+1, b->A2()+1);
        continue;
      }
      BondParmType const& bp = BPA[bpidx];
      double r2 = DIST2_NoImage( fIn.XYZ(b->A1()), fIn.XYZ(b->A2()) );
      double r = sqrt(r2);
      double rdiff = r - bp.Req();
      double ene = bp.Rk() * (rdiff * rdiff);
      Ebond += ene;
#     ifdef DEBUG_ENERGY
      mprintf("\tBond %4u %4i -- %4i: k= %12.5f  x0= %12.5f  r= %12.5f  E= %12.5e\n",
              b - Bonds.begin(), b->A1()+1, b->A2()+1, bp.Rk(), bp.Req(), r, ene);
#     endif
    }
  }
  return Ebond;
}

// -----------------------------------------------------------------------------
/** Angle energy */
double Energy_Amber::E_angle(Frame const& fIn, Topology const& tIn, CharMask const& mask)
{
  // Heavy atom angles
  double Eang = CalcAngleEnergy(fIn, tIn.Angles(), tIn.AngleParm(), mask);
  Eang += CalcAngleEnergy(fIn, tIn.AnglesH(), tIn.AngleParm(), mask);
  return Eang;
}

// Energy_Amber::CalcAngleEnergy()
double Energy_Amber::CalcAngleEnergy(Frame const& fIn, AngleArray const& Angles,
                                     AngleParmArray const& APA, CharMask const& mask)
{
  double Eangle = 0.0;
  for (AngleArray::const_iterator a = Angles.begin(); a != Angles.end(); ++a)
  {
    if (mask.AtomInCharMask(a->A1()) &&
        mask.AtomInCharMask(a->A2()) &&
        mask.AtomInCharMask(a->A3()))
    {
      int apidx = a->Idx();
      if (apidx < 0) {
        if (debug_ > 0)
          mprintf("Warning: Angle %i -- %i -- %i has no parameters.\n",
                   a->A1()+1, a->A2()+1, a->A3()+1);
        continue;
      }
      AngleParmType const& ap = APA[apidx];
      double theta = CalcAngle(fIn.XYZ(a->A1()), fIn.XYZ(a->A2()), fIn.XYZ(a->A3()));
      double tdiff = theta - ap.Teq();
      double ene = ap.Tk() * (tdiff * tdiff);
      Eangle += ene;
#     ifdef DEBUG_ENERGY
      mprintf("\tAngle %4u %4i -- %4i -- %4i: k= %12.5f  x0= %12.5f  t= %12.5f  E= %12.5e\n",
              a - Angles.begin(), a->A1()+1, a->A2()+1, a->A3()+1, ap.Tk(), ap.Teq(), theta, ene);
#     endif
    }
  }
  return Eangle;
}

// -----------------------------------------------------------------------------
/** Dihedral energy */
double Energy_Amber::E_torsion(Frame const& fIn, Topology const& tIn, CharMask const& mask)
{
  // Heavy atom dihedrals
  double Edih = CalcTorsionEnergy(fIn, tIn.Dihedrals(), tIn.DihedralParm(), mask);
  Edih += CalcTorsionEnergy(fIn, tIn.DihedralsH(), tIn.DihedralParm(), mask);
  return Edih;
}

// Energy_Amber::CalcTorsionEnergy()
double Energy_Amber::CalcTorsionEnergy(Frame const& fIn, DihedralArray const& Dihedrals,
                                       DihedralParmArray const& DPA, CharMask const& mask)
{
  double Edih = 0.0;
  for (DihedralArray::const_iterator d = Dihedrals.begin(); d != Dihedrals.end(); d++)
  {
    if (mask.AtomInCharMask(d->A1()) &&
        mask.AtomInCharMask(d->A2()) &&
        mask.AtomInCharMask(d->A3()) &&
        mask.AtomInCharMask(d->A4()))
    {
      int dpidx = d->Idx();
      if (dpidx < 0) {
        if (debug_ > 0)
          mprintf("Warning: Dihedral %i -- %i -- %i -- %i has no parameters.\n",
                   d->A1()+1, d->A2()+1, d->A3()+1, d->A4()+1);
        continue;
      }
      DihedralParmType const& dp = DPA[dpidx];
      double phi = Torsion(fIn.XYZ(d->A1()), fIn.XYZ(d->A2()),
                           fIn.XYZ(d->A3()), fIn.XYZ(d->A4()));
      double ene = dp.Pk() * (1.0 + cos(dp.Pn() * phi - dp.Phase()));
      Edih += ene;
#     ifdef DEBUG_ENERGY
      mprintf("\tDihedral %4u %4i -- %4i -- %4i -- %4i: pk= %12.5f  "
                "pn= %12.5f  phase= %12.5f  p= %12.5f  E= %12.5e\n",
              d - Dihedrals.begin(), d->A1()+1, d->A2()+1, d->A3()+1, d->A4()+1,
              dp.Pk(), dp.Pn(), dp.Phase(), phi, ene);
#     endif
    }
  }
  return Edih;
}

// -----------------------------------------------------------------------------
/** 1-4 nonbond energy */
double Energy_Amber::E_14_Nonbond(Frame const& fIn, Topology const& tIn, CharMask const& mask,
                                  double& Eq14)
{
  Eq14 = 0.0;
  // Heavy atom dihedrals
  double Evdw14 = Calc_14_Energy(fIn, tIn.Dihedrals(), tIn.DihedralParm(), tIn, mask, Eq14);
  Evdw14 += Calc_14_Energy(fIn, tIn.DihedralsH(), tIn.DihedralParm(), tIn, mask, Eq14);
  return Evdw14;
}

// Energy_Amber::Calc_14_Energy()
double Energy_Amber::Calc_14_Energy(Frame const& fIn, DihedralArray const& Dihedrals,
                                    DihedralParmArray const& DPA, Topology const& tIn,
                                    CharMask const& mask, double& Eq14)
{
  double Evdw14 = 0.0;
  for (DihedralArray::const_iterator d = Dihedrals.begin(); d != Dihedrals.end(); d++)
  {
    if (d->Type() == DihedralType::NORMAL &&
        mask.AtomInCharMask(d->A1()) && mask.AtomInCharMask(d->A4()))
    {
      int dpidx = d->Idx();
      if (dpidx < 0) {
        if (debug_ > 0)
          mprintf("Warning: 1-4 pair %i -- %i has no parameters.\n", d->A1()+1, d->A4()+1);
        continue;
      }
      DihedralParmType const& dp = DPA[dpidx];
      double rij2 = DIST2_NoImage( fIn.XYZ(d->A1()), fIn.XYZ(d->A4()) );
      double rij = sqrt( rij2 );
      // VDW
      NonbondType const& LJ = tIn.GetLJparam(d->A1(), d->A4());
      double r2    = 1.0 / rij2;
      double r6    = r2 * r2 * r2;
      double r12   = r6 * r6;
      double f12   = LJ.A() * r12;  // A/r^12
      double f6    = LJ.B() * r6;   // B/r^6
      double e_vdw = f12 - f6;      // (A/r^12)-(B/r^6)
      e_vdw /= dp.SCNB();
      Evdw14 += e_vdw;
      // Coulomb
      double qiqj = QFAC * tIn[d->A1()].Charge() * tIn[d->A4()].Charge();
      double e_elec = qiqj / rij;
      e_elec /= dp.SCEE();
      Eq14 += e_elec;
#     ifdef DEBUG_ENERGY
      mprintf("\tEVDW14  %4i -- %4i: A=  %12.5e  B=  %12.5e  r2= %12.5f  E= %12.5e\n",
              d->A1()+1, d->A4()+1, LJ.A(), LJ.B(), rij2, e_vdw);
      mprintf("\tEELEC14 %4i -- %4i: q1= %12.5e  q2= %12.5e  r=  %12.5f  E= %12.5e\n",
              d->A1()+1, d->A4()+1, tIn[d->A1()].Charge(), tIn[d->A4()].Charge(),
              rij, e_elec);
#     endif
    }
  }
  return Evdw14;
}

// -----------------------------------------------------------------------------
// Energy_Amber::E_Nonbond()
double Energy_Amber::E_Nonbond(Frame const& fIn, Topology const& tIn, AtomMask const& mask,
                               double& EelecOut, ExclusionArray const& Excluded)
{
  double Evdw = 0.0;
  double Eelec = 0.0;
  int idx1;
# ifdef _OPENMP
# pragma omp parallel private(idx1) reduction(+ : Eelec, Evdw)
  {
# pragma omp for
# endif
  for (idx1 = 0; idx1 < mask.Nselected(); idx1++)
  {
    int atom1 = mask[idx1];
    // Set up coord for this atom
    const double* crd1 = fIn.XYZ( atom1 );
    // Set up exclusion list for this atom
    // TODO refactor inner loop to be more like StructureCheck, more efficient.
    ExclusionArray::ExListType::const_iterator excluded_idx = Excluded[idx1].begin();
    for (int idx2 = idx1 + 1; idx2 < mask.Nselected(); idx2++)
    {
      int atom2 = mask[idx2];
      // Advance excluded list up to current selected atom
      while (excluded_idx != Excluded[idx1].end() && *excluded_idx < idx2) ++excluded_idx;
      // If atom is excluded, just increment to next excluded atom.
      if (excluded_idx != Excluded[idx1].end() && idx2 == *excluded_idx)
        ++excluded_idx;
      else {
        double rij2 = DIST2_NoImage( crd1, fIn.XYZ( atom2 ) );
        double rij = sqrt( rij2 );
        // VDW
        NonbondType const& LJ = tIn.GetLJparam(atom1, atom2);
        double r2    = 1.0 / rij2;
        double r6    = r2 * r2 * r2;
        double r12   = r6 * r6;
        double f12   = LJ.A() * r12;  // A/r^12
        double f6    = LJ.B() * r6;   // B/r^6
        double e_vdw = f12 - f6;      // (A/r^12)-(B/r^6)
        Evdw += e_vdw;
        // Coulomb
        double qiqj = QFAC * tIn[atom1].Charge() * tIn[atom2].Charge();
        double e_elec = qiqj / rij;
        Eelec += e_elec;
#       ifdef DEBUG_ENERGY
        mprintf("\tEVDW  %4i -- %4i: A=  %12.5e  B=  %12.5e  r2= %12.5f  E= %12.5e\n",
                atom1+1, atom2+1, LJ.A(), LJ.B(), rij2, e_vdw);
        mprintf("\tEELEC %4i -- %4i: q1= %12.5e  q2= %12.5e  r=  %12.5f  E= %12.5e\n",
                atom1+1, atom2+1, tIn[atom1].Charge(), tIn[atom2].Charge(),
                rij, e_elec);
#       endif
      }
    }
  }
# ifdef _OPENMP
  } // END omp parallel
# endif
  EelecOut = Eelec;
  return Evdw;
}

// -----------------------------------------------------------------------------
double Energy_Amber::E_VDW(Frame const& fIn, Topology const& tIn, AtomMask const& mask,
                           ExclusionArray const& Excluded)
{
  double Evdw = 0.0;
  int idx1;
# ifdef _OPENMP
# pragma omp parallel private(idx1) reduction(+ : Evdw)
  {
# pragma omp for
# endif
  for (idx1 = 0; idx1 < mask.Nselected(); idx1++)
  {
    int atom1 = mask[idx1];
    // Set up coord for this atom
    const double* crd1 = fIn.XYZ( atom1 );
    // Set up exclusion list for this atom
    // TODO refactor inner loop to be more like StructureCheck
    ExclusionArray::ExListType::const_iterator excluded_idx = Excluded[idx1].begin();
    for (int idx2 = idx1 + 1; idx2 < mask.Nselected(); idx2++)
    {
      int atom2 = mask[idx2];
      // Advance excluded list up to current selected atom
      while (excluded_idx != Excluded[idx1].end() && *excluded_idx < idx2) ++excluded_idx;
      // If atom is excluded, just increment to next excluded atom.
      if (excluded_idx != Excluded[idx1].end() && idx2 == *excluded_idx)
        ++excluded_idx;
      else {
        double rij2 = DIST2_NoImage( crd1, fIn.XYZ( atom2 ) );
        // VDW
        NonbondType const& LJ = tIn.GetLJparam(atom1, atom2);
        double r2    = 1.0 / rij2;
        double r6    = r2 * r2 * r2;
        double r12   = r6 * r6;
        double f12   = LJ.A() * r12;  // A/r^12
        double f6    = LJ.B() * r6;   // B/r^6
        double e_vdw = f12 - f6;      // (A/r^12)-(B/r^6)
        Evdw += e_vdw;
#       ifdef DEBUG_ENERGY
        mprintf("\tEVDW  %4i -- %4i: A=  %12.5e  B=  %12.5e  r2= %12.5f  E= %12.5e\n",
                atom1+1, atom2+1, LJ.A(), LJ.B(), rij2, e_vdw);
#       endif
      }
    }
  }
# ifdef _OPENMP
  } // END omp parallel
# endif
  return Evdw;
}

// -----------------------------------------------------------------------------
double Energy_Amber::E_Elec(Frame const& fIn, Topology const& tIn, AtomMask const& mask,
                            ExclusionArray const& Excluded)
{
  double Eelec = 0.0;
  int idx1;
# ifdef _OPENMP
# pragma omp parallel private(idx1) reduction(+ : Eelec)
  {
# pragma omp for
# endif
  for (idx1 = 0; idx1 < mask.Nselected(); idx1++)
  {
    int atom1 = mask[idx1];
    // Set up coord for this atom
    const double* crd1 = fIn.XYZ( atom1 );
    // Set up exclusion list for this atom
    // TODO refactor inner loop to be like StructureCheck
    ExclusionArray::ExListType::const_iterator excluded_idx = Excluded[idx1].begin();
    for (int idx2 = idx1 + 1; idx2 < mask.Nselected(); idx2++)
    {
      int atom2 = mask[idx2];
      // Advance excluded list up to current selected atom
      while (excluded_idx != Excluded[idx1].end() && *excluded_idx < idx2) ++excluded_idx;
      // If atom is excluded, just increment to next excluded atom.
      if (excluded_idx != Excluded[idx1].end() && idx2 == *excluded_idx)
        ++excluded_idx;
      else {
        double rij2 = DIST2_NoImage( crd1, fIn.XYZ( atom2 ) );
        double rij = sqrt( rij2 );
        // Coulomb
        double qiqj = QFAC * tIn[atom1].Charge() * tIn[atom2].Charge();
        double e_elec = qiqj / rij;
        Eelec += e_elec;
#       ifdef DEBUG_ENERGY
        mprintf("\tEELEC %4i -- %4i: q1= %12.5e  q2= %12.5e  r=  %12.5f  E= %12.5e\n",
                atom1+1, atom2+1, tIn[atom1].Charge(), tIn[atom2].Charge(),
                rij, e_elec);
#       endif
      }
    }
  }
# ifdef _OPENMP
  } // END omp parallel
# endif
  return Eelec;
}
// -----------------------------------------------------------------------------
double Energy_Amber::E_DirectSum(Frame const& fIn, Topology const& tIn, AtomMask const& mask,
                                 ExclusionArray const& Excluded,
                                 int n_points)
{
  // Direct sum
  double Edirect = E_Elec(fIn, tIn, mask, Excluded);
  // Sum over images.
  double Eimage = 0.0;
  // Cache npoints values, excluding this cell (0,0,0)
  std::vector<Vec3> Cells;
  int Ncells = (2*n_points)+1;
  Cells.reserve( (Ncells*Ncells*Ncells) - 1 );
  for (int ix = -n_points; ix <= n_points; ix++)
    for (int iy = -n_points; iy <= n_points; iy++)
      for (int iz = -n_points; iz <= n_points; iz++)
        if (ix != 0 || iy != 0 || iz != 0)
          Cells.push_back( Vec3(ix, iy, iz) );
  // Outer loop over atoms (i)
  for (AtomMask::const_iterator atom1 = mask.begin(); atom1 != mask.end(); ++atom1)
  {
//    mprintf("\nDEBUG: Atom %i\n", *atom1+1);
    Vec3 T1( fIn.XYZ(*atom1) );
    // Inner loop over atoms (j)
    for (AtomMask::const_iterator atom2 = mask.begin(); atom2 != mask.end(); ++atom2)
    {
      Vec3 frac2 = fIn.BoxCrd().FracCell() * Vec3(fIn.XYZ( *atom2 )); // atom j in fractional coords
      double qiqj = QFAC * tIn[*atom1].Charge() * tIn[*atom2].Charge();
      // Loop over images of atom j
      for (std::vector<Vec3>::const_iterator ixyz = Cells.begin(); ixyz != Cells.end(); ++ixyz)
      {
//        mprintf("DEBUG: Atom %4i to %4i Image %3i %3i %3i", *atom1+1, *atom2+1, ix, iy, iz);
        // atom j image back in Cartesian space minus atom i in Cartesian space.
        Vec3 dxyz = fIn.BoxCrd().UnitCell().TransposeMult(frac2 + *ixyz) - T1;
        double rij2 = dxyz.Magnitude2();
        double rij = sqrt(rij2);
//        mprintf(" Distance= %g\n", rij);
        double e_elec = qiqj / rij;
        Eimage += e_elec;
      }
    } // atom j
  } // atom i
  return Edirect + (Eimage/2.0);
}

// -----------------------------------------------------------------------------
/** Calculate the kinetic energy using given velocities. */
double Energy_Amber::E_Kinetic(Frame const& fIn, AtomMask const& mask)
{
  if (!fIn.HasVelocity()) return 0.0;
  double ke = 0.0;
  for (AtomMask::const_iterator at = mask.begin(); at != mask.end(); ++at)
  {
    const double* vxyz = fIn.VelXYZ( *at );
    //mprintf("DEBUG: Atom %i vxyz = %12.4f %12.4f %12.4f\n", *at+1, vxyz[0], vxyz[1], vxyz[2]);
    double v2 = vxyz[0]*vxyz[0] + vxyz[1]*vxyz[1] + vxyz[2]*vxyz[2];
    ke += (fIn.Mass(*at) * v2);
  }
  return 0.5 * ke;
}

/** Calculate the kinetic energy assuming a Velocity Verlet / leapfrog scheme,
  * i.e. the given velocities are one half step ahead of the given coordinates
  * and forces.
  * \param fIn Current frame with coordinates, +0.5*dt velocities, and forces.
  * \param mask Selected atoms.
  * \param dt Time step in ps.
  */
double Energy_Amber::E_Kinetic_VV(Frame const& fIn, AtomMask const& mask, double dt)
{
  // Vn = Vh - 0.5 * dt * Fn/m
  // V = velocity, F = force, n = step, h = n + half step, dt = time step, m = mass
  if (!fIn.HasVelocity() || !fIn.HasForce()) return 0.0;
  double ke = 0.0;
  double dthalf = dt * Constants::AMBERTIME_TO_PS * 0.5;
  Vec3 velN(0.0);
  for (AtomMask::const_iterator at = mask.begin(); at != mask.end(); ++at)
  {
    const double* vxyz = fIn.VelXYZ( *at );
    double mass = fIn.Mass( *at );
    //const double* fxyz = fIn.FrcXYZ( *at );
    Vec3 fxyz = Vec3(fIn.FrcXYZ(*at));
    velN[0] = vxyz[0] - dthalf * fxyz[0] / mass; 
    velN[1] = vxyz[1] - dthalf * fxyz[1] / mass; 
    velN[2] = vxyz[2] - dthalf * fxyz[2] / mass; 
    // DEBUG
    /*
    if (at == mask.begin()) {
      double vold[3];
      for (int i = 0; i < 3; i++) {
        vold[i] = vxyz[i] - (dt*Constants::AMBERTIME_TO_PS) * fxyz[i] / tIn[*at].Mass();
        mprintf("DBG: v vold m %6i%12.7f%12.7f%12.7f\n", i+1, vxyz[i], vold[i], tIn[*at].Mass());
      }
      mprintf("\n");
    }
    */
    ke += (mass * velN.Magnitude2());
    //for (int i = 0; i != 3; i++) {
    //  double v2 = (vxyz[i] + vold[i]);
    //  ke += (tIn[*at].Mass() * 0.25 * (v2 * v2));
    //}
  }
  return 0.5 * ke;
}

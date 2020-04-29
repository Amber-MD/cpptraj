#include "PotentialTerm_LJ_Coulomb.h"
#include "Topology.h"
#include "CharMask.h"

/**  CONSTRUCTOR */
PotentialTerm_LJ_Coulomb::PotentialTerm_LJ_Coulomb() :
  PotentialTerm(SIMPLE_LJ_Q),
  nonbond_(0),
  E_vdw_(0),
  E_elec_(0)
{}

/** Setup nonbonds. */
int PotentialTerm_LJ_Coulomb::SetupTerm(Topology const& topIn, CharMask const& maskIn,
                                  EnergyArray& Earray)
{
  selectedAtoms_.clear();
  nonSelectedAtoms_.clear();

  selectedAtoms_.reserve( maskIn.Nselected() );
  nonSelectedAtoms_.reserve( topIn.Natom() - maskIn.Nselected() );
  for (int i = 0; i < topIn.Natom(); i++) {
    if (maskIn.AtomInCharMask( i ))
      selectedAtoms_.push_back( i );
    else
      nonSelectedAtoms_.push_back( i );
  }

  atoms_ = &(topIn.Atoms());
  nonbond_ = &(topIn.Nonbond());

  return 0;
}

/** Calculate LJ and coulomb forces between specified atoms. */
static inline void NonBondKernel(Frame& frameIn, int idx, int jdx,
                                 double LJA, double LJB,
                                 CharMask const& maskIn,
                                 double& E_vdw, double& E_elec)
{
  const double* XYZ0 = frameIn.XYZ( idx );
  const double* XYZ1 = frameIn.XYZ( jdx );
  double rx = XYZ0[0] - XYZ1[0];
  double ry = XYZ0[1] - XYZ1[1];
  double rz = XYZ0[2] - XYZ1[2];
  double rij2 = rx*rx + ry*ry + rz*rz;
  if (rij2 > 0) {
    //double rij;
    //if (rij2 < nbcut2) {
    //  rij2 = nbcut2;
    //  // Make rij really big to scale down the coulomb part.
    //  rij = 99999;
    //} else
    //  rij = sqrt( rij2 );
    double rij = sqrt( rij2 );
    //double dfx = 0;
    //double dfy = 0;
    //double dfz = 0;
    // VDW
    double r2    = 1.0 / rij2;
    double r6    = r2 * r2 * r2;
    double r12   = r6 * r6;
    double f12   = LJA * r12;  // A/r^12
    double f6    = LJB * r6;   // B/r^6
    double e_vdw = f12 - f6;   // (A/r^12)-(B/r^6)
    E_vdw += e_vdw;
    // VDW force
    double fvdw = ((12*f12) - (6*f6)) * r2; // (12A/r^13)-(6B/r^7)
    double dfx = rx * fvdw;
    double dfy = ry * fvdw;
    double dfz = rz * fvdw;
    // COULOMB
    double qiqj = 1; // Give each atom charge of 1
    double e_elec = 1.0 * (qiqj / rij); // 1.0 is electrostatic constant, not really needed
    E_elec += e_elec;
    // COULOMB force
    double felec = e_elec / rij; // kes * (qiqj / r) * (1/r)
    dfx += rx * felec;
    dfy += ry * felec;
    dfz += rz * felec;
    // Apply forces
    if (maskIn.AtomInCharMask(idx)) {
      double* fxyz = frameIn.fAddress() + (3*idx);
      fxyz[0] += dfx;
      fxyz[1] += dfy;
      fxyz[2] += dfz;
    }
    if (maskIn.AtomInCharMask(jdx)) {
      double* fxyz = frameIn.fAddress() + (3*jdx);
      fxyz[0] -= dfx;
      fxyz[1] -= dfy;
      fxyz[2] -= dfz;
    }
  } // END rij > 0
}

static inline NonbondType const& GetLJparam(std::vector<Atom> const& atoms,
                                            NonbondParmType const& nonbond, int a1, int a2)
{
  int nbindex = nonbond.GetLJindex( atoms[a1].TypeIndex(), atoms[a2].TypeIndex() );
  //if (nbindex < 0) // Means Amber Hbond, return A = B = 0.0
  //  return LJ_EMPTY;
  return nonbond.NBarray( nbindex );
}

/** Calculate nonbonds. */
void PotentialTerm_LJ_Coulomb::CalcForce(Frame& frameIn, CharMask const& maskIn) const {
  *E_vdw_ = 0.0;
  *E_elec_ = 0.0;
  // First loop is each non-selected atom to each selected atom.
  // There is no overlap between the two.
  for (Iarray::const_iterator idx = nonSelectedAtoms_.begin();
                              idx != nonSelectedAtoms_.end(); ++idx)
  {
    for (Iarray::const_iterator jdx = selectedAtoms_.begin(); jdx != selectedAtoms_.end(); ++jdx)
    {
      // Ignore if idx and jdx are bonded.
      if (!(*atoms_)[*idx].IsBondedTo(*jdx))
      {
        NonbondType LJ = GetLJparam(*atoms_, *nonbond_, *idx, *jdx);
        NonBondKernel(frameIn, *idx, *jdx, LJ.A(), LJ.B(), maskIn, *E_vdw_, *E_elec_);
      } // END idx not bonded to jdx
    } // END inner loop over jdx
  } // END outer loop over idx
  // Second loop is each selected atom to each other selected atom.
  for (Iarray::const_iterator idx0 = selectedAtoms_.begin();
                              idx0 != selectedAtoms_.end(); ++idx0)
  {
    for (Iarray::const_iterator idx1 = idx0 + 1; idx1 != selectedAtoms_.end(); ++idx1)
    {
      // Ignore if idx0 and idx1 are bonded.
      if (!(*atoms_)[*idx0].IsBondedTo(*idx1))
      {
        NonbondType LJ = GetLJparam(*atoms_, *nonbond_, *idx0, *idx1);
        NonBondKernel(frameIn, *idx0, *idx1, LJ.A(), LJ.B(), maskIn, *E_vdw_, *E_elec_);
      }
    } // END inner loop over idx1
  } // END outer loop over idx0
}

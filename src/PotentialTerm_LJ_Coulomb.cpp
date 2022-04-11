#include "PotentialTerm_LJ_Coulomb.h"
#include "Topology.h"
#include "CharMask.h"
#include "EnergyArray.h"
#include "EnergyKernel_NonBond_Simple.h"
#include "MdOpts.h"
#include "CpptrajStdio.h"
#include "Constants.h" // SMALL
#include <cmath> //fabs

/**  CONSTRUCTOR */
PotentialTerm_LJ_Coulomb::PotentialTerm_LJ_Coulomb() :
  PotentialTerm(SIMPLE_LJ_Q),
  nonbond_(0),
  E_vdw_(0),
  E_elec_(0),
  QFAC_(0),
  cutoff2_(0),
  nExclude_(4)
{}

/** Init nonbonds. */
int PotentialTerm_LJ_Coulomb::InitTerm(MdOpts const& optsIn) {
  nExclude_ = optsIn.N_Exclude();
  QFAC_ = optsIn.CoulombFactor();
  cutoff2_ = optsIn.CutEE();
  if ( fabs(cutoff2_ - optsIn.CutNB()) > Constants::SMALL ) {
    mprinterr("Error: Simple nonbond term does not yet support separate elec/LJ cutoff.\n");
    return 1;
  }
  // Store cutoff^2
  cutoff2_ *= cutoff2_;
  return 0;
}

/** Setup nonbonds. */
int PotentialTerm_LJ_Coulomb::SetupTerm(Topology const& topIn,Box const& boxIn, CharMask const& maskIn,
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
  // Set up exclusion array to ignore self and have full list
  // since the non-selected to selected loop may not have atoms in order.
  AtomMask fullSystem(0, topIn.Natom());
  if (Excluded_.SetupExcluded(topIn.Atoms(), fullSystem, nExclude_,
                              ExclusionArray::EXCLUDE_SELF,
                              ExclusionArray::FULL))
  {
    mprinterr("Error: PotentialTerm_LJ_Coulomb: Could not set up exclusion list.\n");
    return 1;
  }

  atoms_ = &(topIn.Atoms());
  nonbond_ = &(topIn.Nonbond());
  E_vdw_ = Earray.AddType(EnergyArray::E_VDW);
  E_elec_ = Earray.AddType(EnergyArray::E_COULOMB);

  return 0;
}

/** Get LJ parameter for interaction between two atoms. */
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
  EnergyKernel_NonBond_Simple<double> nonbond(cutoff2_);
  // First loop is each non-selected atom to each selected atom.
  // There is no overlap between the two.
  for (Iarray::const_iterator idx = nonSelectedAtoms_.begin();
                              idx != nonSelectedAtoms_.end(); ++idx)
  {
    // Exclusion list for this atom
    ExclusionArray::ExListType const& excluded = Excluded_[*idx];
    for (Iarray::const_iterator jdx = selectedAtoms_.begin(); jdx != selectedAtoms_.end(); ++jdx)
    {
      // Check exclusion
      if (excluded.find( *jdx ) == excluded.end())
      {
        NonbondType LJ = GetLJparam(*atoms_, *nonbond_, *idx, *jdx);
        //NonBondKernel(frameIn, *idx, *jdx, LJ.A(), LJ.B(), maskIn, *E_vdw_, *E_elec_);
        nonbond.Calc_F_E(frameIn, *idx, *jdx, LJ.A(), LJ.B(), QFAC_, (*atoms_)[*idx].Charge(), (*atoms_)[*jdx].Charge(), 1, 1, maskIn, *E_vdw_, *E_elec_);
      } // END idx not bonded to jdx
    } // END inner loop over jdx
  } // END outer loop over idx
  // Second loop is each selected atom to each other selected atom.
  for (Iarray::const_iterator idx0 = selectedAtoms_.begin();
                              idx0 != selectedAtoms_.end(); ++idx0)
  {
    // Exclusion list for this atom
    ExclusionArray::ExListType const& excluded = Excluded_[*idx0];
    for (Iarray::const_iterator idx1 = idx0 + 1; idx1 != selectedAtoms_.end(); ++idx1)
    {
      // Check exclusion
      if (excluded.find( *idx1 ) == excluded.end())
      {
        NonbondType LJ = GetLJparam(*atoms_, *nonbond_, *idx0, *idx1);
        nonbond.Calc_F_E(frameIn, *idx0, *idx1, LJ.A(), LJ.B(), QFAC_, (*atoms_)[*idx0].Charge(), (*atoms_)[*idx1].Charge(), 1, 1, maskIn, *E_vdw_, *E_elec_);
        //NonBondKernel(frameIn, *idx0, *idx1, LJ.A(), LJ.B(), maskIn, *E_vdw_, *E_elec_);
      }
    } // END inner loop over idx1
  } // END outer loop over idx0
}

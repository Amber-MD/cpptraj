#include "PotentialTerm_LJ_Coulomb.h"
#include "Topology.h"
#include "CharMask.h"
#include "EnergyArray.h"
#include "EnergyKernel_NonBond_Simple.h"

/**  CONSTRUCTOR */
PotentialTerm_LJ_Coulomb::PotentialTerm_LJ_Coulomb() :
  PotentialTerm(SIMPLE_LJ_Q),
  nonbond_(0),
  E_vdw_(0),
  E_elec_(0)
{}

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
  // FIXME not actually using charges here.
  EnergyKernel_NonBond_Simple<double> nonbond;
  // First loop is each non-selected atom to each selected atom.
  // There is no overlap between the two.
  for (Iarray::const_iterator idx = nonSelectedAtoms_.begin();
                              idx != nonSelectedAtoms_.end(); ++idx)
  {
    for (Iarray::const_iterator jdx = selectedAtoms_.begin(); jdx != selectedAtoms_.end(); ++jdx)
    {
      // Ignore if idx and jdx are bonded. TODO use exclusion
      if (!(*atoms_)[*idx].IsBondedTo(*jdx))
      {
        NonbondType LJ = GetLJparam(*atoms_, *nonbond_, *idx, *jdx);
        //NonBondKernel(frameIn, *idx, *jdx, LJ.A(), LJ.B(), maskIn, *E_vdw_, *E_elec_);
        nonbond.Calc_F_E(frameIn, *idx, *jdx, LJ.A(), LJ.B(), 1, 1, 1, 1, 1, maskIn, *E_vdw_, *E_elec_);
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
        nonbond.Calc_F_E(frameIn, *idx0, *idx1, LJ.A(), LJ.B(), 1, 1, 1, 1, 1, maskIn, *E_vdw_, *E_elec_);
        //NonBondKernel(frameIn, *idx0, *idx1, LJ.A(), LJ.B(), maskIn, *E_vdw_, *E_elec_);
      }
    } // END inner loop over idx1
  } // END outer loop over idx0
}

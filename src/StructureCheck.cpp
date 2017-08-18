#include <cmath> // sqrt
#include "StructureCheck.h"
#include "CpptrajStdio.h"
#ifdef _OPENMP
#  include <omp.h>
#endif

/// CONSTRUCTOR
StructureCheck::StructureCheck() :
  bondoffset_(1.15),
  nonbondcut2_(0.64), // 0.8^2
  // NOTE: Default of 4.0 Ang for cutoff is from trial and error; seems
  //       to give a good balance between speed and grid size.
  plcut_(4.0),
  checkType_(NO_PL_1_MASK),
  bondcheck_(true),
  saveProblems_(false)
{}

// StructureCheck::SetOptions()
int StructureCheck::SetOptions(bool imageOn, bool checkBonds, bool saveProblemsIn,
                               std::string const& mask1, std::string const& mask2,
                               double overlapCut, double bondLengthOffset, double pairListCut)
{
  image_.InitImaging( imageOn );
  bondcheck_ = checkBonds;
  saveProblems_ = saveProblemsIn;
  bondoffset_ = bondLengthOffset;
  nonbondcut2_ = overlapCut * overlapCut; // Save cutoff squared.
  plcut_ = pairListCut;
  Mask1_.SetMaskString( mask1 );
  if (!mask2.empty())
    Mask2_.SetMaskString( mask2 );
  // TODO error checking
  # ifdef _OPENMP
  // Each thread needs space to store problems to prevent clashes
# pragma omp parallel
  {
# pragma omp master
  {
  thread_problemAtoms_.resize( omp_get_num_threads() );
  }
  }
# endif
  return 0;
}

/** Set up bond arrays in a sorted list for easy access during loop
  * over all pairs of atoms. Only use bonds for which both atoms are in
  * the mask.
  */
void StructureCheck::ProcessBondArray(BondArray const& Bonds, BondParmArray const& Parm,
                                      CharMask const& cMask)
{
  for (BondArray::const_iterator bnd = Bonds.begin(); bnd != Bonds.end(); ++bnd)
  {
    if ( cMask.AtomInCharMask(bnd->A1()) && cMask.AtomInCharMask(bnd->A2()) ) {
      if (bnd->Idx() < 0)
        mprintf("Warning: Bond parameters not present for atoms %i-%i, skipping.\n",
                bnd->A1()+1, bnd->A2()+1);
      else {
        double Req_off = Parm[ bnd->Idx() ].Req() + bondoffset_;
        bondList_.push_back( Problem(bnd->A1(), bnd->A2(), Req_off*Req_off) );
      }
    }
  }
}

/** Set up bond parameters for bonds for which both atoms present in mask. */
void StructureCheck::SetupBondList(AtomMask const& iMask, Topology const& top) {
  CharMask cMask( iMask.ConvertToCharMask(), iMask.Nselected() );
 
  ProcessBondArray(top.Bonds(),  top.BondParm(), cMask);
  ProcessBondArray(top.BondsH(), top.BondParm(), cMask);
}

// StructureCheck::Setup()
int StructureCheck::Setup(Topology const& topIn, Box const& boxIn)
{
  image_.SetupImaging( boxIn.Type() );
  bondList_.clear();
  // Set up first mask
  if ( topIn.SetupIntegerMask( Mask1_ ) ) return 1;
  Mask1_.MaskInfo();
  if (Mask1_.None()) {
    mprinterr("Error: Mask '%s' has no atoms.\n", Mask1_.MaskString());
    return 1;
  }
  checkType_ = NO_PL_1_MASK;
  // Set up bonds if specified.
  if (bondcheck_) SetupBondList(Mask1_, topIn);
  // Set up second mask if specified.
  if ( Mask2_.MaskStringSet() ) {
    if (topIn.SetupIntegerMask( Mask2_ ) ) return 1;
    Mask2_.MaskInfo();
    if (Mask2_.None()) {
      mprinterr("Error: Mask '%s' has no atoms.\n", Mask2_.MaskString());
      return 1;
    }
    int common = Mask1_.NumAtomsInCommon( Mask2_ );
    if (common > 0)
      mprintf("Warning: '%s' has %i atoms in common with '%s'. Some problems may be reported\n"
              "Warning:   more than once.\n", Mask1_.MaskString(), common, Mask2_.MaskString());
    // Outer mask should be the one with the most atoms.
    if ( Mask2_.Nselected() > Mask1_.Nselected() ) {
      OuterMask_ = Mask2_;
      InnerMask_ = Mask1_;
    } else {
      OuterMask_ = Mask1_;
      InnerMask_ = Mask2_;
    }
    if (bondcheck_) SetupBondList(Mask2_, topIn);
    checkType_ = NO_PL_2_MASKS;
  }
  // Check if pairlist should be used.
  if (image_.ImagingEnabled() && !Mask2_.MaskStringSet()) {
    if (pairList_.InitPairList( plcut_, 0.1, 0 )) return 1;
    Matrix_3x3 ucell, recip;
    boxIn.ToRecip(ucell, recip);
    if (pairList_.SetupPairList( boxIn.Type(), boxIn.RecipLengths(recip) )) return 1;
    mprintf("\tUsing pair list.\n");
    checkType_ = PL_1_MASK;
  }
  return 0;
}

// StructureCheck::CheckBonds()
int StructureCheck::CheckBonds(Frame const& currentFrame)
{
  int Nproblems = 0;
  problemAtoms_.clear();
  int bond_max = (int)bondList_.size();
  int idx;
# ifdef _OPENMP
  int mythread;
# pragma omp parallel private(idx,mythread) reduction(+: Nproblems)
  {
  mythread = omp_get_thread_num();
  thread_problemAtoms_[mythread].clear();
# pragma omp for
# endif
  for (idx = 0; idx < bond_max; idx++) {
    double D2 = DIST2_NoImage( currentFrame.XYZ(bondList_[idx].A1()),
                               currentFrame.XYZ(bondList_[idx].A2()) );
    if (D2 > bondList_[idx].D()) {
      ++Nproblems;
      if (saveProblems_) {
#       ifdef _OPENMP
        thread_problemAtoms_[mythread]
#       else
        problemAtoms_
#       endif
          .push_back(Problem(bondList_[idx].A1(), bondList_[idx].A2(), sqrt(D2)));
      }
    }
  } // END loop over bonds
# ifdef _OPENMP
  } // END pragma omp parallel
# endif

  return Nproblems;   
}

#ifdef _OPENMP
/** Combine problems from each thread into problemAtoms_ */
void StructureCheck::ConsolidateProblems() {
  for (unsigned int thread = 0; thread != thread_problemAtoms_.size(); ++thread)
    for (Parray::const_iterator p = thread_problemAtoms_[thread].begin();
                                p != thread_problemAtoms_[thread].end(); ++p)
      problemAtoms_.push_back( *p );
}
#endif

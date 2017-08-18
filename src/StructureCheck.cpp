#include <cmath> // sqrt
#include <algorithm> // sort
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
  // Sort bond list
  if (bondcheck_) std::sort(bondList_.begin(), bondList_.end());
  return 0;
}

/** Sort problem list. If OpenMP, combine problems from each thread
  * into problemAtoms_ first.
  */
void StructureCheck::ConsolidateProblems() {
# ifdef _OPENMP
  for (unsigned int thread = 0; thread != thread_problemAtoms_.size(); ++thread)
    for (Parray::const_iterator p = thread_problemAtoms_[thread].begin();
                                p != thread_problemAtoms_[thread].end(); ++p)
      problemAtoms_.push_back( *p );
# endif
  std::sort( problemAtoms_.begin(), problemAtoms_.end() );
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
  ConsolidateProblems();

  return Nproblems;   //TODO problemAtoms_.size()?
}

/** Check for bad overlaps; use a pair list to speed things up.
  * \return Number of bad overlaps.
  */
int StructureCheck::PL1_CheckOverlap(Frame const& currentFrame, Matrix_3x3 const& ucell,
                                     Matrix_3x3 const& recip)
{
  int Nproblems = 0;
  pairList_.CreatePairList(currentFrame, ucell, recip, Mask1_);
  problemAtoms_.clear();

  int cidx;
# ifdef _OPENMP
  int mythread;
# pragma omp parallel private(cidx,mythread) reduction(+: Nproblems)
  {
  mythread = omp_get_thread_num();
  thread_problemAtoms_[mythread].clear();
# pragma omp for
# endif 
  for (cidx = 0; cidx < pairList_.NGridMax(); cidx++)
  {
    PairList::CellType const& thisCell = pairList_.Cell( cidx );
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
        // Calc interaction of atom to all other atoms in thisCell.
        for (PairList::CellType::const_iterator it1 = it0 + 1;
                                                it1 != thisCell.end(); ++it1)
        {
          Vec3 const& xyz1 = it1->ImageCoords();
          Vec3 dxyz = xyz1 - xyz0;
          double D2 = dxyz.Magnitude2();
          if (D2 < nonbondcut2_) {
            ++Nproblems;
            if (saveProblems_) {
#             ifdef _OPENMP
              thread_problemAtoms_[mythread]
#             else
              problemAtoms_
#             endif
                .push_back(Problem(Mask1_[it0->Idx()], Mask1_[it1->Idx()], sqrt(D2)));
            }
          }
        } // END loop over all other atoms in thisCell
        // Loop over all neighbor cells
        for (unsigned int nidx = 1; nidx != cellList.size(); nidx++)
        {
          PairList::CellType const& nbrCell = pairList_.Cell( cellList[nidx] );
          // Translate vector for neighbor cell
          Vec3 const& tVec = pairList_.TransVec( transList[nidx] );
          // Loop over every atom in nbrCell
          for (PairList::CellType::const_iterator it1 = nbrCell.begin();
                                                  it1 != nbrCell.end(); ++it1)
          {
            Vec3 const& xyz1 = it1->ImageCoords();
            Vec3 dxyz = xyz1 + tVec - xyz0;
            double D2 = dxyz.Magnitude2();
            if (D2 < nonbondcut2_) {
              ++Nproblems;
              if (saveProblems_) {
#               ifdef _OPENMP
                thread_problemAtoms_[mythread]
#               else
                problemAtoms_
#               endif
                  .push_back(Problem(Mask1_[it0->Idx()], Mask1_[it1->Idx()], sqrt(D2)));
              }
            }
          } // END loop over atoms in neighbor cell
        } // END loop over neighbor cells
      } // END loop over atoms in thisCell
    } // END cell not empty
  } // END loop over cells
# ifdef _OPENMP
  } // END omp parallel
# endif
  ConsolidateProblems();

  return Nproblems;
}

// StructureCheck::Mask2_CheckOverlap()
int StructureCheck::Mask2_CheckOverlap(Frame const& currentFrame, Matrix_3x3 const& ucell,
                                       Matrix_3x3 const& recip)
{
  int Nproblems = 0;
  // Calculation of all atoms in Mask1 to all atoms in Mask2
  int outer_max = OuterMask_.Nselected();
  int inner_max = InnerMask_.Nselected();
  int nmask1;
# ifdef _OPENMP
  int mythread;
# pragma omp parallel private(mythread,nmask1) reduction(+: Nproblems)
  {
  mythread = omp_get_thread_num();
  thread_problemAtoms_[mythread].clear();
# pragma omp for
# endif
  for (nmask1 = 0; nmask1 < outer_max; nmask1++) {
    int atom1 = OuterMask_[nmask1];
    for (int nmask2 = 0; nmask2 < inner_max; nmask2++) {
      int atom2 = InnerMask_[nmask2];
      if (atom1 != atom2) {
        double D2 = DIST2( currentFrame.XYZ(atom1), currentFrame.XYZ(atom2),
                           image_.ImageType(), currentFrame.BoxCrd(), ucell, recip);
        if (D2 < nonbondcut2_) {
          ++Nproblems;
          if (saveProblems_) {
#           ifdef _OPENMP
            thread_problemAtoms_[mythread]
#           else
            problemAtoms_
#           endif
              .push_back(Problem(atom1, atom2, sqrt(D2)));
          }
        }
      }
    } // END loop over inner mask
  } // END loop over outer mask
# ifdef _OPENMP
  } // END pragma omp parallel
# endif
  ConsolidateProblems();

  return Nproblems;
}

// StructureCheck::Mask1_CheckOverlap()
int StructureCheck::Mask1_CheckOverlap(Frame const& currentFrame, Matrix_3x3 const& ucell,
                                       Matrix_3x3 const& recip)
{
  int Nproblems = 0;
  // Calculation of atoms in Mask1 to all other atoms in Mask1
  int mask1_max = Mask1_.Nselected();
  int nmask1;
# ifdef _OPENMP
  int mythread;
# pragma omp parallel private(mythread,nmask1) reduction(+: Nproblems)
  {
  mythread = omp_get_thread_num();
  thread_problemAtoms_[mythread].clear();
# pragma omp for schedule(dynamic)
# endif
  for (nmask1 = 0; nmask1 < mask1_max; nmask1++) {
    int atom1 = Mask1_[nmask1];
    for (int nmask2 = nmask1 + 1; nmask2 < mask1_max; nmask2++) {
      int atom2 = Mask1_[nmask2];
      double D2 = DIST2( currentFrame.XYZ(atom1), currentFrame.XYZ(atom2),
                         image_.ImageType(), currentFrame.BoxCrd(), ucell, recip);
      if (D2 < nonbondcut2_) {
        ++Nproblems;
        if (saveProblems_) {
#           ifdef _OPENMP
            thread_problemAtoms_[mythread]
#           else
            problemAtoms_
#           endif
              .push_back(Problem(atom1, atom2, sqrt(D2)));
        }
      }
    } // END inner loop over Mask1
  } // END outer loop over Mask1
# ifdef _OPENMP
  } // END pragma omp parallel
# endif
  ConsolidateProblems();

  return Nproblems;
}

//  StructureCheck::CheckOverlaps()
int StructureCheck::CheckOverlaps(Frame const& currentFrame) {
  Matrix_3x3 ucell, recip;
  if (checkType_ == PL_1_MASK || image_.ImageType() == NONORTHO)
    currentFrame.BoxCrd().ToRecip(ucell, recip);
  int Nproblems = 0;
  switch (checkType_) {
    case PL_1_MASK     : Nproblems = PL1_CheckOverlap(currentFrame, ucell, recip); break;
    case NO_PL_2_MASKS : Nproblems = Mask2_CheckOverlap(currentFrame, ucell, recip); break;
    case NO_PL_1_MASK  : Nproblems = Mask1_CheckOverlap(currentFrame, ucell, recip); break;
  }
  return Nproblems;
}

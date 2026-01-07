#include <cmath> // sqrt
#include <algorithm> // sort
#include "StructureCheck.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "CharMask.h"
#include "DistRoutines.h"
#include "Structure/LeastSquaresPlane.h"
#include "Topology.h"
#ifdef _OPENMP
#  include <omp.h>
#endif

/// CONSTRUCTOR
StructureCheck::StructureCheck() :
  XP_Exclude_Mask_("&!@/XP"),
  bondoffset_(1.15),
  bondMinOffset_(0.50),
  nonbondcut2_(0.64), // (0.8 Ang)^2
  plcut_(8.0),
  ring_shortd2_(0.25), // Short ring-bond 0.5 Ang. distance cutoff
  ring_dcut2_(1.562500), // 1.25 Ang. ring-bond distance cutoff
  ring_acut_(1.1344640138), // 65 deg. ringvec-bondvec angle cutoff
  checkType_(NO_PL_1_MASK),
  debug_(0),
  bondcheck_(true),
  ringcheck_(true),
  saveProblems_(false),
  checkExtraPts_(false),
  lastFmt_(F_ATOM)
{}

/** \return Ring-bond short distance cutoff in Ang. */
double StructureCheck::RingShortDist() const {
  return sqrt(ring_shortd2_);
}

/** \return Ring-bond distance cutoff in Ang. */
double StructureCheck::RingDist() const {
  return sqrt(ring_dcut2_);
}

/** \return Ring vector-bond vector angle cutoff in degrees. */
double StructureCheck::RingAngleCut_Deg() const {
  return ring_acut_ * Constants::RADDEG;
}

/** Create mask string based on if we are excluding extra points or not. */
std::string StructureCheck::checkMaskStr(std::string const& maskIn)
const
{
  std::string maskOut;
  if (maskIn.empty())
    maskOut = "*" + XP_Exclude_Mask_;
  else
    maskOut = maskIn + XP_Exclude_Mask_;
  mprintf("DEBUG: StructureCheck::checkMaskStr(): String is %s\n", maskOut.c_str());
  return maskOut;
}

/** Set options with default for extra points (exclude) */
int StructureCheck::SetOptions(bool imageOn, bool checkBonds, bool saveProblemsIn,
                               int debugIn,
                               std::string const& mask1, std::string const& mask2, 
                               double overlapCut, double bondLengthOffset, double bondMinOffset,
                               double pairListCut,
                               bool checkRings, double ring_shortd, double ring_dcut, double ring_acut)
{
  return SetOptions(imageOn, checkBonds, saveProblemsIn, false, debugIn, mask1, mask2, "", overlapCut,
                    bondLengthOffset, bondMinOffset, pairListCut, checkRings, ring_shortd, ring_dcut, ring_acut);
}

// StructureCheck::SetOptions()
int StructureCheck::SetOptions(bool imageOn, bool checkBonds, bool saveProblemsIn, bool checkExtraPtsIn,
                               int debugIn,
                               std::string const& mask1, std::string const& mask2, std::string const& XP_Exclude_MaskIn,
                               double overlapCut, double bondLengthOffset, double bondMinOffset,
                               double pairListCut,
                               bool checkRings, double ring_shortd, double ring_dcut, double ring_acut)
{
  imageOpt_.InitImaging( imageOn );
  bondcheck_ = checkBonds;
  ringcheck_ = checkRings;
  saveProblems_ = saveProblemsIn;
  checkExtraPts_ = checkExtraPtsIn;
  debug_ = debugIn;
  bondoffset_ = bondLengthOffset;
  bondMinOffset_ = bondMinOffset;
  nonbondcut2_ = overlapCut * overlapCut; // Save cutoff squared.
  plcut_ = pairListCut;
  if (ring_shortd > 0)
    ring_shortd2_ = ring_shortd * ring_shortd;
  if (ring_dcut > 0)
    ring_dcut2_ = ring_dcut * ring_dcut;
  if (ring_acut > 0)
    ring_acut_ = ring_acut * Constants::DEGRAD;
  if (!XP_Exclude_MaskIn.empty())
    XP_Exclude_Mask_ = XP_Exclude_MaskIn;
  else
    XP_Exclude_Mask_ = "&!@/XP"; // Select extra points by element
  if (checkExtraPts_) {
    XP_Exclude_Mask_.clear();
  }
  if (Mask1_.SetMaskString( checkMaskStr(mask1) )) return 1;
  if (!mask2.empty()) {
    if (Mask2_.SetMaskString( checkMaskStr(mask2) )) return 1;
  }
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
  rings_.SetDebug( debug_ );
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
        double Req_off = Parm[ bnd->Idx() ].Req();// + bondoffset_;
        //bondList_.push_back( Problem(bnd->A1(), bnd->A2(), Req_off*Req_off) );
        bondList_.push_back( Problem(bnd->A1(), bnd->A2(), Req_off) );
      }
    }
  }
}

/** Set up bond parameters for bonds for which both atoms present in mask. */
void StructureCheck::SetupBondList(AtomMask const& iMask, Topology const& top) {
  CharMask cMask( iMask.ConvertToCharMask(), iMask.Nselected() );
 
  ProcessBondArray(top.Bonds(),  top.BondParm(), cMask);
  ProcessBondArray(top.BondsH(), top.BondParm(), cMask);
  // Only look at bonds to heavy atoms for rings
  ringBonds_.clear();
  ringBonds_.reserve( top.Bonds().size() );
  for (BondArray::const_iterator bnd = top.Bonds().begin(); bnd != top.Bonds().end(); ++bnd)
  {
    if ( cMask.AtomInCharMask(bnd->A1()) && cMask.AtomInCharMask(bnd->A2()) ) {
      ringBonds_.push_back( Btype(bnd->A1(), bnd->A2()) );
    }
  }
}

// StructureCheck::Setup()
int StructureCheck::Setup(Topology const& topIn, Box const& boxIn)
{
  imageOpt_.SetupImaging( boxIn.HasBox() );
  bondList_.clear();
  if (!checkExtraPts_)
    mprintf("\tExcluding extra points using mask expression '%s'.\n", XP_Exclude_Mask_.c_str());
  // Set up first mask
  if ( topIn.SetupIntegerMask( Mask1_ ) ) return 1;
  if (Mask1_.None()) {
    mprinterr("Error: Mask '%s' has no atoms.\n", Mask1_.MaskString());
    return 1;
  }
  checkType_ = NO_PL_1_MASK;
  // Set up bonds if specified.
  if (bondcheck_ || ringcheck_) SetupBondList(Mask1_, topIn);
  // Set up second mask if specified.
  if ( Mask2_.MaskStringSet() ) {
    if (topIn.SetupIntegerMask( Mask2_ ) ) return 1;
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
    if (bondcheck_ || ringcheck_) SetupBondList(Mask2_, topIn);
    checkType_ = NO_PL_2_MASKS;
  }
  // Check if pairlist should be used.
  ExclusionArray::SelfOpt ex_self_opt = ExclusionArray::NO_EXCLUDE_SELF;
  ExclusionArray::ListOpt ex_list_opt = ExclusionArray::ONLY_GREATER_IDX;
  if (imageOpt_.ImagingEnabled() && !Mask2_.MaskStringSet()) {
    mprintf("\tUsing pair list.\n");
    // If pairlist cutoff < 0 try to use a heuristic to estimate.
    if (plcut_ < 0) {
      // Get the minimum dimension length
      double minDimLen = boxIn.Param(Box::X);
      if (boxIn.Param(Box::Y) < minDimLen)
        minDimLen = boxIn.Param(Box::Y);
      if (boxIn.Param(Box::Z) < minDimLen)
        minDimLen = boxIn.Param(Box::Z);
      // Have cutoff be sqrt of min dim len
      plcut_ = sqrt(minDimLen);
      if (plcut_ < 8.0) {
        plcut_ = 8.0;
        mprintf("\tSet pairlist cutoff to default of %f Ang.\n", plcut_);
      } else
        mprintf("\tSet pairlist cutoff to %f Ang. based on minimum box dimension %f Ang.\n",
                plcut_, minDimLen);
    }
    if (pairList_.InitPairList( plcut_, 0.1, debug_ )) {
      mprinterr("Error: StructureCheck: Could not init pair list.\n");
      return 1;
    }
    if (pairList_.SetupPairList( boxIn )) {
      mprinterr("Error: StructureCheck: Could not setup pair list.\n");
      return 1;
    }
    checkType_ = PL_1_MASK;
    ex_self_opt = ExclusionArray::EXCLUDE_SELF;
    ex_list_opt = ExclusionArray::FULL;
  }
  // Set up exclusion list
  if (checkType_ == PL_1_MASK || checkType_ == NO_PL_1_MASK) {
    mprintf("\tExcluding bond interactions.\n");
    // Set up atom exclusion list. Distance of 2 since we are not
    // yet checking angles and dihedrals. 
    if (Excluded_.SetupExcluded(topIn.Atoms(), Mask1_, 2, ex_self_opt, ex_list_opt)) {
      mprinterr("Error: StructureCheck: Could not set up excluded atoms list.\n");
      return 1;
    }
  } else
    mprintf("\tNot using exclusions.\n");
  // Sort bond list
  if (bondcheck_ || ringcheck_) std::sort(bondList_.begin(), bondList_.end());
  // Find rings
  if (ringcheck_) {
    if (rings_.SetupRingFinder(topIn, Mask1_)) {
      mprinterr("Error: Could not set up ring finder.\n");
      return 1;
    }
    if (debug_ > 0) rings_.PrintRings(topIn);
  }

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
    double dist = sqrt(D2);
    //if (D2 > bondList_[idx].D()) {
    if (dist > (bondList_[idx].D()+bondoffset_) ||
        dist < (bondList_[idx].D()-bondMinOffset_))
    {
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
  lastFmt_ = F_BOND;

  return Nproblems;
}

/** Check if any bonds are passing through rings. */
int StructureCheck::CheckRings(Frame const& currentFrame) {
  return CheckRings(currentFrame, rings_, ringBonds_);
}

/// DEBUG - For printing ring atoms in mask
//static inline void printRingAtoms(AtomMask const& maskIn) {
//  for (AtomMask::const_iterator it = maskIn.begin(); it != maskIn.end(); ++it)
//    mprintf(" %i", *it);
//}

/** Check ring and bond */
void StructureCheck::ring_bond_check(int& Nproblems,
                                     double dist2, Btype const& bnd, Vec3 const& vbond,
                                     AtomMask const& ringMask,
                                     Cpptraj::Structure::LeastSquaresPlane const& ringVec
#                                    ifdef _OPENMP
                                     , int mythread
#                                    endif
                                    )
{
  bool ring_intersect = false;
  // Bond intersects ring if it meets the short cutoff or if the angle
  // between the bond and the ring normal is less than a cutoff.
  if (dist2 < ring_shortd2_) {
    //mprintf("DEBUG: Bond %i - %i near ring (%f).",
    //        bnd.A1(), bnd.A2(), sqrt(dist2));
    //printRingAtoms(ringMask);
    //mprintf("\n");
    ring_intersect = true;
  } else {
    //mprintf("DEBUG: Bond %i - %i near ring (%f), doing angle check.",
    //        bnd.A1(), bnd.A2(), sqrt(dist2));
    //printRingAtoms(ringMask);
    //mprintf("\n");
    // Get the angle
    double ang_in_rad = vbond.Angle( ringVec.Nxyz() );
    // Wrap the angle between 0-90 degrees
    if (ang_in_rad > Constants::PIOVER2)
      ang_in_rad = Constants::PI - ang_in_rad;
    //mprintf("DEBUG:\t\tWrapped angle is %f deg.\n", Constants::RADDEG*ang_in_rad);
    if (ang_in_rad < ring_acut_) {
      //mprintf("DEBUG: Bond %i - %i near ring (%f) Ang= %f deg.",
      //        bnd.A1(), bnd.A2(), sqrt(dist2),
      //        Constants::RADDEG*ang_in_rad);
      ////printRingAtoms(ringMask);
      //mprintf("\n");
      ring_intersect = true;
    }
  }
  if (ring_intersect) {
    //mprintf("DEBUG: Bond intersects ring.\n");
    ++Nproblems;
    if (saveProblems_) {
      // Do not use constructor since we do not want to sort atoms
      Problem newProb;
      newProb.SetProb( bnd.A1(), ringMask.back(), sqrt(dist2) );
#     ifdef _OPENMP
      thread_problemAtoms_[mythread]
#     else
      problemAtoms_
#     endif
        .push_back( newProb );
    }
  } // END angle cutoff satisfied
}

bool StructureCheck::check_bond_not_in_ring(AtomMask const& ringMask, Btype const& bnd)
{
  return ( !ringMask.IsSelected(bnd.A1()) &&
           !ringMask.IsSelected(bnd.A2()) );
}

/** Check if any bonds are passing through rings, use pairlist. */
int StructureCheck::checkRings_PL(Frame const& currentFrame,
                                  Cpptraj::Structure::RingFinder const& rings,
                                  std::vector<Btype> const& ringBonds,
                                  std::vector<Cpptraj::Structure::LeastSquaresPlane> const& RingVecs)
{
  int Nproblems = 0;
  // Need to create a pseudo-frame containing bond midpoints and
  // ring centers.
  unsigned int pseudo_natom = ringBonds.size() + rings.Nrings();

  Frame pseudoFrame( pseudo_natom );
  pseudoFrame.SetBox( currentFrame.BoxCrd() );
  //pseudoFrame.ClearAtoms();

  // Bonds first
  int idx = 0;
  int bond_max = (int)ringBonds.size();
  std::vector<Vec3> Vbonds( ringBonds.size() );
# ifdef _OPENMP
# pragma omp parallel private(idx)
  {
# pragma omp for
# endif
  for (idx = 0; idx < bond_max; idx++)
  {
    const double* xyz1 = currentFrame.XYZ( ringBonds[idx].A1() );
    const double* xyz2 = currentFrame.XYZ( ringBonds[idx].A2() );
    Vec3 vbond( xyz2[0] - xyz1[0],
                xyz2[1] - xyz1[1],
                xyz2[2] - xyz1[2] );
    Vec3 vmid = Vec3(xyz1) + (vbond * 0.5);
    vbond.Normalize();
    pseudoFrame.SetXYZ(idx, vmid);
    Vbonds[idx] = vbond;
  }
# ifdef _OPENMP
  } // END pragma omp parallel
# endif

  // Rings second
  idx = bond_max;
  for (unsigned int jdx = 0; jdx != rings.Nrings(); jdx++, idx++) {
    Cpptraj::Structure::LeastSquaresPlane const& ringVec = RingVecs[jdx];
    pseudoFrame.SetXYZ(idx, ringVec.Cxyz());
  }

  // Pseudo mask selecting everything
  AtomMask pseudoMask(0, pseudo_natom);

  // Set up the pair list
  int retVal = pairList_.CreatePairList(pseudoFrame,
                                        pseudoFrame.BoxCrd().UnitCell(),
                                        pseudoFrame.BoxCrd().FracCell(), pseudoMask);
  if (retVal < 0) {
    // Treat grid setup failure as one problem.
    mprinterr("Error: Grid setup for ring-bond check failed.\n");
    return 1;
  } else if (retVal > 0) {
    // Atoms off the grid should count as problems as well.
    Nproblems = retVal;
  }

  // Loop over pairlist cells
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
        bool it0_is_bond = (it0->Idx() < bond_max);
        // Calc interaction of atom to all other atoms in thisCell.
        for (PairList::CellType::const_iterator it1 = it0 + 1;
                                                it1 != thisCell.end(); ++it1)
        {
          bool it1_is_bond = (it1->Idx() < bond_max);
          if (it0_is_bond != it1_is_bond) {
            AtomMask const* ringMask = 0;
            Cpptraj::Structure::LeastSquaresPlane const* ringVec = 0;
            Btype const* rbnd = 0;
            Vec3 const* vbond = 0;
            if (it0_is_bond) {
              ringMask = &(rings[it1->Idx() - bond_max]);
              ringVec = &(RingVecs[it1->Idx() - bond_max]);
              rbnd = &(ringBonds[it0->Idx()]);
              vbond = &(Vbonds[it0->Idx()]);
            } else {
              ringMask = &(rings[it0->Idx() - bond_max]);
              ringVec = &(RingVecs[it0->Idx() - bond_max]);
              rbnd = &(ringBonds[it1->Idx()]);
              vbond = &(Vbonds[it1->Idx()]);
            }
            if (check_bond_not_in_ring( *ringMask, *rbnd )) {
              Vec3 const& xyz1 = it1->ImageCoords();
              Vec3 dxyz = xyz1 - xyz0;
              double dist2 = dxyz.Magnitude2();
              if (dist2 < ring_dcut2_) {
                ring_bond_check(Nproblems, dist2, *rbnd, *vbond, *ringMask, *ringVec
#                               ifdef _OPENMP
                                , mythread
#                               endif
                               );
              }  // END outer distance cutoff satisfied
            } // END bond not in ring
          } // END it0 and it1 are not the same type
        } // END loop over all other atoms of thisCell.
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
            bool it1_is_bond = (it1->Idx() < bond_max);
            if (it0_is_bond != it1_is_bond) {
              AtomMask const* ringMask = 0;
              Cpptraj::Structure::LeastSquaresPlane const* ringVec = 0;
              Btype const* rbnd = 0;
              Vec3 const* vbond = 0;
              if (it0_is_bond) {
                ringMask = &(rings[it1->Idx() - bond_max]);
                ringVec = &(RingVecs[it1->Idx() - bond_max]);
                rbnd = &(ringBonds[it0->Idx()]);
                vbond = &(Vbonds[it0->Idx()]);
              } else {
                ringMask = &(rings[it0->Idx() - bond_max]);
                ringVec = &(RingVecs[it0->Idx() - bond_max]);
                rbnd = &(ringBonds[it1->Idx()]);
                vbond = &(Vbonds[it1->Idx()]);
              }
              if (check_bond_not_in_ring( *ringMask, *rbnd )) {
                Vec3 const& xyz1 = it1->ImageCoords();
                Vec3 dxyz = xyz1 + tVec - xyz0;
                double dist2 = dxyz.Magnitude2();
                if (dist2 < ring_dcut2_) {
                  ring_bond_check(Nproblems, dist2, *rbnd, *vbond, *ringMask, *ringVec
#                                 ifdef _OPENMP
                                  , mythread
#                                 endif
                                 );
                } // END outer distance cutoff satisfied
              } // END bond not in ring
            } // END it0 and it1 are not the same type
          } // END loop over every atom in neighbor cell
        } // END loop over all neighbor cells
      } // END loop over all atoms of thisCell
    } // END thisCell has atoms
  } // END loop over pairlist cells
# ifdef _OPENMP
  } // END omp parallel
# endif

  return Nproblems;
}

/** Check if any bonds are passing through rings, no pairlist. */
int StructureCheck::checkRings_NoPL(Frame const& currentFrame,
                                    Cpptraj::Structure::RingFinder const& rings,
                                    std::vector<Btype> const& ringBonds,
                                    std::vector<Cpptraj::Structure::LeastSquaresPlane> const& RingVecs)
{
  int Nproblems = 0;
  int idx = 0;
  int ring_max = (int)rings.Nrings();
  // Loop over bonds
  int bond_max = (int)ringBonds.size();
# ifdef _OPENMP
  int mythread;
# pragma omp parallel private(idx,mythread) reduction(+: Nproblems)
  {
  mythread = omp_get_thread_num();
  thread_problemAtoms_[mythread].clear();
# pragma omp for
# endif
  for (idx = 0; idx < bond_max; idx++)
  {
    const double* xyz1 = currentFrame.XYZ( ringBonds[idx].A1() );
    const double* xyz2 = currentFrame.XYZ( ringBonds[idx].A2() );
    Vec3 vbond( xyz2[0] - xyz1[0],
                xyz2[1] - xyz1[1],
                xyz2[2] - xyz1[2] );
    Vec3 vmid = Vec3(xyz1) + (vbond * 0.5);
    vbond.Normalize();
    // Loop over rings
    for (int jdx = 0; jdx < ring_max; jdx++)
    {
      // Make sure this bond is not in this ring
      AtomMask const& ringMask = rings[jdx];
      if ( !ringMask.IsSelected(ringBonds[idx].A1()) &&
           !ringMask.IsSelected(ringBonds[idx].A2()) )
      {
        // Get the center distance
        Cpptraj::Structure::LeastSquaresPlane const& ringVec = RingVecs[jdx];
        double dist2 = DIST2_NoImage( vmid.Dptr(), ringVec.Cxyz().Dptr() );
        if (dist2 < ring_dcut2_) {
          ring_bond_check(Nproblems, dist2, ringBonds[idx], vbond, ringMask, ringVec
#                         ifdef _OPENMP
                          , mythread
#                         endif
                         );
/*
          bool ring_intersect = false;
          // Bond intersects ring if it meets the short cutoff or if the angle
          // between the bond and the ring normal is less than a cutoff.
          if (dist2 < ring_shortd2_) {
            mprintf("DEBUG: Bond %i - %i near ring %i (%f).",
                    ringBonds[idx].A1(), ringBonds[idx].A2(), jdx, sqrt(dist2));
            printRingAtoms(ringMask);
            mprintf("\n");
            ring_intersect = true;
          } else {
            mprintf("DEBUG: Bond %i - %i near ring %i (%f), doing angle check.",
                    ringBonds[idx].A1(), ringBonds[idx].A2(), jdx, sqrt(dist2));
            printRingAtoms(ringMask);
            mprintf("\n");
            // Get the angle
            double ang_in_rad = vbond.Angle( ringVec.Nxyz() );
            // Wrap the angle between 0-90 degrees
            if (ang_in_rad > Constants::PIOVER2)
              ang_in_rad = Constants::PI - ang_in_rad;
            mprintf("DEBUG:\t\tWrapped angle is %f deg.\n", Constants::RADDEG*ang_in_rad);
            if (ang_in_rad < ring_acut_) {
              mprintf("DEBUG: Bond %i - %i near ring %i (%f) Ang= %f deg.",
                      ringBonds[idx].A1(), ringBonds[idx].A2(), jdx, sqrt(dist2),
                      Constants::RADDEG*ang_in_rad);
              //printRingAtoms(ringMask);
              mprintf("\n");
              ring_intersect = true;
            }
          }
          if (ring_intersect) {
            mprintf("DEBUG: Bond intersects ring.\n");
            ++Nproblems;
            if (saveProblems_) {
              // Do not use constructor since we do not want to sort atoms
              Problem newProb;
              newProb.SetProb( ringBonds[idx].A1(), ringMask.back(), sqrt(dist2) );
#             ifdef _OPENMP
              thread_problemAtoms_[mythread]
#             else
              problemAtoms_
#             endif
                .push_back( newProb );
            }
          } // END angle cutoff satisfied
*/
        } // END distance cutoff satisfied
      } // END both bond atoms are not part of this ring
    } // END loop over rings
  } // END loop over bonds
# ifdef _OPENMP
  } // END pragma omp parallel
# endif
  return Nproblems;
}

/** Check if any bonds are passing through rings. */
int StructureCheck::CheckRings(Frame const& currentFrame, Cpptraj::Structure::RingFinder const& rings, std::vector<Btype> const& ringBonds)
{
  problemAtoms_.clear();
  // Get all ring vectors
  typedef std::vector<Cpptraj::Structure::LeastSquaresPlane> Larray;
  Larray RingVecs;
  RingVecs.resize( rings.Nrings() );
  int idx = 0;
  int ring_max = (int)rings.Nrings();
# ifdef _OPENMP
# pragma omp parallel private(idx)
  {
# pragma omp for
# endif
  for (idx = 0; idx < ring_max; idx++)
  {
    // false = use geometric center
    RingVecs[idx].CalcLeastSquaresPlane( currentFrame, rings[idx], false );
  }
# ifdef _OPENMP
  } // END pragma omp parallel
# endif

  int Nproblems = 0;
  if (checkType_ == PL_1_MASK) {
    // TODO check for valid box
    Nproblems = checkRings_PL(currentFrame, rings,  ringBonds, RingVecs);
  } else {
    Nproblems = checkRings_NoPL(currentFrame, rings,  ringBonds, RingVecs);
  }

  ConsolidateProblems();
  lastFmt_ = F_RING;

  return Nproblems;
}

/** Check for bad overlaps; use a pair list to speed things up.
  * \return Number of bad overlaps.
  */
int StructureCheck::PL1_CheckOverlap(Frame const& currentFrame)
{
  int Nproblems = 0;
  int retVal = pairList_.CreatePairList(currentFrame,
                                        currentFrame.BoxCrd().UnitCell(),
                                        currentFrame.BoxCrd().FracCell(), Mask1_);
  if (retVal < 0) {
    // Treat grid setup failure as one problem.
    mprinterr("Error: Grid setup failed.\n");
    return 1;
  } else if (retVal > 0) {
    // Atoms off the grid should count as problems as well.
    Nproblems = retVal;
  }
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
        // Exclusion list for this atom
        ExclusionArray::ExListType const& excluded = Excluded_[it0->Idx()];
        // Calc interaction of atom to all other atoms in thisCell.
        for (PairList::CellType::const_iterator it1 = it0 + 1;
                                                it1 != thisCell.end(); ++it1)
        {
          // If atom not excluded, calculate distance
          if (excluded.find(it1->Idx()) == excluded.end())
          {
            Vec3 const& xyz1 = it1->ImageCoords();
            Vec3 dxyz = xyz1 - xyz0;
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
          } // END atom not excluded
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
            // If atom not excluded, calculate distance
            if (excluded.find(it1->Idx()) == excluded.end())
            {
              Vec3 const& xyz1 = it1->ImageCoords();
              Vec3 dxyz = xyz1 + tVec - xyz0;
              double D2 = dxyz.Magnitude2();
              if (D2 < nonbondcut2_) {
                ++Nproblems;
                if (saveProblems_) {
#                 ifdef _OPENMP
                  thread_problemAtoms_[mythread]
#                 else
                  problemAtoms_
#                 endif
                    .push_back(Problem(Mask1_[it0->Idx()], Mask1_[it1->Idx()], sqrt(D2)));
                }
              }
            } // END atom not excluded
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

/** Check for and record non-bonded clashes. */
void StructureCheck::DistanceCheck(Frame const& currentFrame, int atom1, int atom2,
                                   Parray& problemAtoms, int& Nproblems)
const
{
  double D2 = DIST2( imageOpt_.ImagingType(), currentFrame.XYZ(atom1), currentFrame.XYZ(atom2),
                     currentFrame.BoxCrd());
  if (D2 < nonbondcut2_) {
    ++Nproblems;
    if (saveProblems_) {
      problemAtoms.push_back(Problem(atom1, atom2, sqrt(D2)));
    }
  }
}

// StructureCheck::Mask2_CheckOverlap()
int StructureCheck::Mask2_CheckOverlap(Frame const& currentFrame)
{
  problemAtoms_.clear();
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
        DistanceCheck(currentFrame, atom1, atom2,
#                     ifdef _OPENMP
                      thread_problemAtoms_[mythread],
#                     else
                      problemAtoms_,
#                     endif
                      Nproblems);
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
int StructureCheck::Mask1_CheckOverlap(Frame const& currentFrame)
{
  problemAtoms_.clear();
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
    // Inner loop is broken into 2 phases. First the atoms excluded from
    // interacting with atom1 are handled, then everything else. This
    // leverages the fact that excluded atoms tend to be close in
    // sequence.
    int nmask2 = nmask1 + 1;
    int atom2 = Mask1_[nmask2];
    // Advance excluded list up to current selected atom
    ExclusionArray::ExListType::const_iterator ex = Excluded_[nmask1].begin();
    while (ex != Excluded_[nmask1].end() && *ex < nmask2) ++ex;
    // Continue looping over excluded atoms until there are no more to exclude.
    while (ex != Excluded_[nmask1].end()) {
      atom2 = Mask1_[nmask2];
      if (nmask2 == *ex)
        ++ex;
      else {
        DistanceCheck(currentFrame, atom1, atom2,
#                     ifdef _OPENMP
                      thread_problemAtoms_[mythread],
#                     else
                      problemAtoms_,
#                     endif
                      Nproblems);
      }
      nmask2++;
    } // End loop over atom1 exclusion list
    // Now, no more interactions to exclude.
    for (; nmask2 < mask1_max; nmask2++) {
      atom2 = Mask1_[nmask2];
      DistanceCheck(currentFrame, atom1, atom2,
#                   ifdef _OPENMP
                    thread_problemAtoms_[mythread],
#                   else
                    problemAtoms_,
#                   endif
                    Nproblems);
    }
  } // END outer loop over Mask1
# ifdef _OPENMP
  } // END pragma omp parallel
# endif
  ConsolidateProblems();

  return Nproblems;
}

//  StructureCheck::CheckOverlaps()
int StructureCheck::CheckOverlaps(Frame const& currentFrame) {
  //mprintf("Entering CheckOverlaps()\n");
  int Nproblems = 0;

  // First check the box if imaging is enabled.
  bool box_is_bad = false;
  if (imageOpt_.ImagingEnabled()) {
    if (currentFrame.BoxCrd().CheckBox() != Box::BOX_OK) {
      // Box is bad
      box_is_bad = true;
      mprintf("Warning: Disabling imaging due to problem with box.\n");
      Nproblems++;
      // Since there is a problem with the box, disable imaging.
      if (imageOpt_.ImagingEnabled())
        imageOpt_.SetImageType( ImageOption::NO_IMAGE );
    } else {
      // Box is ok
      imageOpt_.SetImageType( currentFrame.BoxCrd().Is_X_Aligned_Ortho() );
    }
  }

  if (checkType_ == PL_1_MASK) {
    // First, check the box before using pairlist.
    if (box_is_bad) {
      // Since there is a problem with the box, do not use the pair list
      mprintf("Warning: Not using pair list due to problem with box.\n");
      Nproblems += Mask1_CheckOverlap(currentFrame);
    } else {
      //mprintf("DEBUG: Box is ok.\n");
      Nproblems += PL1_CheckOverlap(currentFrame);
    }
  } else if (checkType_ == NO_PL_2_MASKS) {
    Nproblems = Mask2_CheckOverlap(currentFrame);
  } else { // (checkType_ == NO_PL_1_MASK)
    Nproblems = Mask1_CheckOverlap(currentFrame);
  }
  //mprintf("Exiting CheckOverlaps() with %i\n\n", Nproblems);
  lastFmt_ = F_ATOM;

  return Nproblems;
}

/** Write formats for problems. */
const char* StructureCheck::Fmt_[] = {
  "%i\t Warning: Atoms %i:%s and %i:%s are close (%.2f)\n",    ///< F_ATOM
  "%i\t Warning: Unusual bond length %i:%s to %i:%s (%.2f)\n", ///< F_BOND
  "%i\t Warning: Bond involving atom %i:%s intersects ring involving atom %i:%s (%.2f)\n" ///< F_RING
};

/** Write current problems to the given file. */
void StructureCheck::WriteProblemsToFile(CpptrajFile* outfile, int frameNum, Topology const& top) const {
  if (outfile == 0) return;
  for (const_iterator p = begin(); p != end(); ++p) {
    outfile->Printf(Fmt_[lastFmt_], frameNum,
                    p->A1()+1, top.TruncResAtomName(p->A1()).c_str(),
                    p->A2()+1, top.TruncResAtomName(p->A2()).c_str(), p->D());
  }
}

#include <cmath> //sqrt
#include <algorithm> // sort
#include "Action_CheckStructure.h"
#include "CpptrajStdio.h"
#ifdef _OPENMP
#  include <omp.h>
#endif

// CONSTRUCTOR
Action_CheckStructure::Action_CheckStructure() :
  bondoffset_(1.15),
  nonbondcut2_(0.64), // 0.8^2
  // NOTE: Default of 4.0 Ang for cutoff is from trial and error; seems
  //       to give a good balance between speed and grid size.
  plcut_(4.0),
  outfile_(0),
  CurrentParm_(0),
  num_problems_(0),
  silent_(false),
  skipBadFrames_(false),
  bondcheck_(true),
  usePairList_(false)
{}

void Action_CheckStructure::Help() const {
  mprintf("\t[<mask>] [around <mask2>] [reportfile <report>] [noimage] [skipbadframes]\n"
          "\t[offset <offset>] [cut <cut>] [nobondcheck] [silent]\n"
          "  Check frames for atomic overlaps and unusual bond lengths\n");
}

// Action_CheckStructure::SeparateInit()
int Action_CheckStructure::SeparateInit(bool imageOn, std::string const& mask1,
                                        std::string const& mask2, std::string const& fname,
                                        double cutIn, double offsetIn, bool silentIn,
                                        DataFileList& DFL)
{
  image_.InitImaging( imageOn );
  bondoffset_ = offsetIn;
  nonbondcut2_ = cutIn * cutIn; // Save cutoff squared.
  silent_ = silentIn;
  if (!silent_)
    outfile_ = DFL.AddCpptrajFile(fname, "Structure check", DataFileList::TEXT, true);
  Mask1_.SetMaskString( mask1 );
  if (!mask2.empty())
    Mask2_.SetMaskString( mask2 );
  return 0;
}

// Action_CheckStructure::Init()
Action::RetType Action_CheckStructure::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get Keywords
  std::string around = actionArgs.GetStringKey("around");
  SeparateInit( !(actionArgs.hasKey("noimage")), actionArgs.GetMaskNext(),
                around, actionArgs.GetStringKey("reportfile"),
                actionArgs.getKeyDouble("cut",0.8),
                actionArgs.getKeyDouble("offset",1.15),
                actionArgs.hasKey("silent"), init.DFL() );
  // DoAction-only keywords.
  plcut_ = actionArgs.getKeyDouble("plcut", 4.0);
  bondcheck_ = !actionArgs.hasKey("nobondcheck");
  skipBadFrames_ = actionArgs.hasKey("skipbadframes");
  DataFile* dfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  num_problems_ = init.DSL().AddSet( DataSet::INTEGER, actionArgs.GetStringNext(), "CHECK" );
  if (num_problems_ == 0) return Action::ERR;
  if (dfile != 0) dfile->AddDataSet( num_problems_ );
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

  mprintf("    CHECKSTRUCTURE: Checking atoms in mask '%s'",Mask1_.MaskString());
  if (Mask2_.MaskStringSet())
    mprintf(" around mask '%s'", Mask2_.MaskString());
  if (!image_.UseImage())
    mprintf(", imaging off");
  if (outfile_ != 0)
    mprintf(", output to %s", outfile_->Filename().full());
  mprintf(".\n");
  mprintf("\tNumber of problems in each frame will be saved to set '%s'\n",
          num_problems_->legend());
  if (dfile != 0)
    mprintf("\tNumber of problems each frame will be written to '%s'\n",
            dfile->DataFilename().full());
  if (!bondcheck_) {
    mprintf("\tChecking inter-atomic distances only.\n");
    mprintf("\tWarnings will be printed for non-bond distances < %.2f Ang.\n", sqrt(nonbondcut2_));
  } else {
    mprintf("\tChecking inter-atomic and bond distances.\n");
    mprintf("\tWarnings will be printed for bond lengths > eq + %.2f Ang\n",
            bondoffset_);
    mprintf("\tand non-bond distances < %.2f Ang.\n", sqrt(nonbondcut2_));
  }
  if (skipBadFrames_) {
    mprintf("\tFrames with problems will be skipped.\n");
#   ifdef MPI
    if (init.TrajComm().Size() > 1)
      mprintf("Warning: Skipping frames in parallel can cause certain actions "
                       "(e.g. 'rms') to hang.\n"
              "Warning:   In addition, trajectories written after skipping "
                       "frames may have issues.\n");
#   endif
  }
  if (silent_)
    mprintf("\tStructure warning messages will be suppressed.\n");
  mprintf("\tCutoff for building pair list is %f Ang.\n", plcut_);
# ifdef _OPENMP
  mprintf("\tParallelizing calculation with %zu threads.\n", thread_problemAtoms_.size());
# endif
  return Action::OK;
}

/** Set up bond arrays in a sorted list for easy access during loop
  * over all pairs of atoms. Only use bonds for which both atoms are in
  * the mask.
  */
void Action_CheckStructure::ProcessBondArray(BondArray const& Bonds, BondParmArray const& Parm,
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
void Action_CheckStructure::SetupBondList(AtomMask const& iMask, Topology const& top) {
  CharMask cMask( iMask.ConvertToCharMask(), iMask.Nselected() );
 
  ProcessBondArray(top.Bonds(),  top.BondParm(), cMask);
  ProcessBondArray(top.BondsH(), top.BondParm(), cMask);
}

// Action_CheckStructure::SeparateSetup()
int Action_CheckStructure::SeparateSetup(Topology const& top, Box::BoxType btype, bool checkBonds)
{
  image_.SetupImaging( btype );
  bondList_.clear();
  // Set up masks
  if ( top.SetupIntegerMask( Mask1_ ) ) return 1;
  Mask1_.MaskInfo();
  if (Mask1_.None()) {
    mprinterr("Error: Mask '%s' has no atoms.\n", Mask1_.MaskString());
    return 1;
  }
  if (checkBonds) SetupBondList(Mask1_, top);
  if ( Mask2_.MaskStringSet() ) {
    if (top.SetupIntegerMask( Mask2_ ) ) return 1;
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
    if (checkBonds) SetupBondList(Mask2_, top);
  }
  return 0;
}

// Action_CheckStructure::Setup()
Action::RetType Action_CheckStructure::Setup(ActionSetup& setup) {
  CurrentParm_ = setup.TopAddress();
  if (SeparateSetup( setup.Top(), setup.CoordInfo().TrajBox().Type(),bondcheck_ ))
    return Action::ERR;
  if (bondcheck_)
    mprintf("\tChecking %u bonds.\n", bondList_.size());
  // Print imaging info for this parm
  if (image_.ImagingEnabled())
    mprintf("\tImaging on.\n");
  else
    mprintf("\timaging off.\n");
  // Set up pairlist. Only use if imaging and not using second mask.
  usePairList_ = false;
  if (image_.ImagingEnabled() && !Mask2_.MaskStringSet()) {
    if (pairList_.InitPairList( plcut_, 0.1, 0 )) return Action::ERR;
    Matrix_3x3 ucell, recip;
    Box const& box = setup.CoordInfo().TrajBox();
    box.ToRecip(ucell, recip);
    if (pairList_.SetupPairList( box.Type(), box.RecipLengths(recip) )) return Action::ERR;
    mprintf("\tUsing pair list.\n");
    usePairList_ = true;
  }
  return Action::OK;
}

/** Check for bad bond lengths. */
int Action_CheckStructure::CheckBonds(int frameNum, Frame const& currentFrame, Topology const& top)
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
      if (outfile_ != 0) {
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
  if (outfile_ != 0) WriteProblems(BondFmt_, frameNum, top); 

  return Nproblems;   
}

/** Check for bad overlaps; use a pair list to speed things up.
  * \return Number of bad overlaps.
  */
int Action_CheckStructure::PL_CheckOverlap(int frameNum, Frame const& currentFrame,
                                           Topology const& top)
{
  int Nproblems = 0;
  Matrix_3x3 ucell, recip; // ToFrac, ToCart
  currentFrame.BoxCrd().ToRecip(ucell, recip);
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
            if (outfile_ != 0) {
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
              if (outfile_ != 0) {
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
  if (outfile_ != 0) WriteProblems(AtomFmt_, frameNum, top); 

  return Nproblems;
}

const char* Action_CheckStructure::BondFmt_ =
  "%i\t Warning: Unusual bond length %i:%s to %i:%s (%.2f)\n";

const char* Action_CheckStructure::AtomFmt_ =
  "%i\t Warning: Atoms %i:%s and %i:%s are close (%.2f)\n";

/** Consolidate problems from different threads if necessary and write out. */
void Action_CheckStructure::WriteProblems(const char* fmt, int frameNum, Topology const& top) {
# ifdef _OPENMP
  for (unsigned int thread = 0; thread != thread_problemAtoms_.size(); ++thread)
    for (Parray::const_iterator p = thread_problemAtoms_[thread].begin();
                                p != thread_problemAtoms_[thread].end(); ++p)
      problemAtoms_.push_back( *p );
# endif
  std::sort( problemAtoms_.begin(), problemAtoms_.end() );
  for (Parray::const_iterator p = problemAtoms_.begin(); p != problemAtoms_.end(); ++p)
    //outfile_->Printf("%i\t Warning: Atoms %i:%s and %i:%s are close (%.2f)\n", frameNum,
    outfile_->Printf(fmt, frameNum,
                    p->A1()+1, top.TruncResAtomName(p->A1()).c_str(),
                    p->A2()+1, top.TruncResAtomName(p->A2()).c_str(), p->D());
}

/** Check for bad overlaps. */
int Action_CheckStructure::CheckOverlap(int frameNum, Frame const& currentFrame, Topology const& top)
{
  double D2;
  Matrix_3x3 ucell, recip; // ToFrac, ToCart
  int nmask1, nmask2;
  int atom1, atom2;
  int Nproblems = 0;
  problemAtoms_.clear();

  // Get imaging info for non-orthogonal box // TODO Check volume
  if (image_.ImageType()==NONORTHO)
    currentFrame.BoxCrd().ToRecip(ucell, recip);
  if ( Mask2_.MaskStringSet() ) {
    // Calculation of all atoms in Mask1 to all atoms in Mask2
    int outer_max = OuterMask_.Nselected();
    int inner_max = InnerMask_.Nselected();
#   ifdef _OPENMP
    int mythread;
#   pragma omp parallel private(mythread,nmask1,nmask2,atom1,atom2,D2) reduction(+: Nproblems)
    {
    mythread = omp_get_thread_num();
    thread_problemAtoms_[mythread].clear();
#   pragma omp for
#   endif
    for (nmask1 = 0; nmask1 < outer_max; nmask1++) {
      atom1 = OuterMask_[nmask1];
      for (nmask2 = 0; nmask2 < inner_max; nmask2++) {
        atom2 = InnerMask_[nmask2];
        if (atom1 != atom2) {
          D2 = DIST2( currentFrame.XYZ(atom1), currentFrame.XYZ(atom2),
                      image_.ImageType(), currentFrame.BoxCrd(), ucell, recip);
          if (D2 < nonbondcut2_) {
            ++Nproblems;
            if (outfile_ != 0) {
#             ifdef _OPENMP
              thread_problemAtoms_[mythread]
#             else
              problemAtoms_
#             endif
                .push_back(Problem(atom1, atom2, sqrt(D2)));
            }
          }
        }
      } // END loop over inner mask
    } // END loop over outer mask
#   ifdef _OPENMP
    } // END pragma omp parallel
#   endif
  } else {
    // Calculation of atoms in Mask1 to all other atoms in Mask1
    int mask1_max = Mask1_.Nselected();
#   ifdef _OPENMP
    int mythread;
#   pragma omp parallel private(mythread,nmask1,nmask2,atom1,atom2,D2) reduction(+: Nproblems)
    {
    mythread = omp_get_thread_num();
    thread_problemAtoms_[mythread].clear();
#   pragma omp for schedule(dynamic)
#   endif
    for (nmask1 = 0; nmask1 < mask1_max; nmask1++) {
      atom1 = Mask1_[nmask1];
      for (nmask2 = nmask1 + 1; nmask2 < mask1_max; nmask2++) {
        atom2 = Mask1_[nmask2];
        D2 = DIST2( currentFrame.XYZ(atom1), currentFrame.XYZ(atom2),
                    image_.ImageType(), currentFrame.BoxCrd(), ucell, recip);
        if (D2 < nonbondcut2_) {
          ++Nproblems;
          if (outfile_ != 0) {
#             ifdef _OPENMP
              thread_problemAtoms_[mythread]
#             else
              problemAtoms_
#             endif
                .push_back(Problem(atom1, atom2, sqrt(D2)));
          }
        }
      } // END inner loop over Mask1
    } // END outer loop over Mask1
#   ifdef _OPENMP
    } // END pragma omp parallel
#   endif
  }
  if (outfile_ != 0) WriteProblems(AtomFmt_, frameNum, top);

  return Nproblems;
}

// Action_CheckStructure::DoAction()
Action::RetType Action_CheckStructure::DoAction(int frameNum, ActionFrame& frm) {
  int total_problems;
  if (usePairList_)
    total_problems = PL_CheckOverlap(frameNum+1, frm.Frm(), *CurrentParm_);
  else
    total_problems = CheckOverlap(frameNum+1, frm.Frm(), *CurrentParm_);
  if (bondcheck_)
    total_problems += CheckBonds(frameNum+1, frm.Frm(), *CurrentParm_);
  num_problems_->Add( frameNum, &total_problems );
  if (total_problems > 0 && skipBadFrames_)
    return Action::SUPPRESS_COORD_OUTPUT;
  return Action::OK;
}

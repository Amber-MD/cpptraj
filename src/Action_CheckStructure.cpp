#include <algorithm> // min and max

#include "Action_CheckStructure.h"
#include "CpptrajStdio.h"
#ifdef MPI
# include "DataSet_integer.h"
# include "DataSet_double.h"
# include "DataSet_string.h"
#endif

// CONSTRUCTOR
Action_CheckStructure::Action_CheckStructure() :
  outfile_(0),
  CurrentParm_(0),
  num_problems_(0),
  silent_(false),
  skipBadFrames_(false)
{}

void Action_CheckStructure::Help() const {
  mprintf("\t[<mask>] [around <mask2>] [reportfile <report>] [noimage]\n"
          "\t[skipbadframes] [offset <offset>] [minoffset <minoffset>]\n"
          "\t[cut <cut>] [nobondcheck] [silent] [plcut <cut>]\n"
          "\t[{noringcheck | [ringshortdist <rsdist>] [ringdcut <ringdcut>]\n"
          "\t                [ringacut <ringacut>]}]\n"
          "\t[{checkxp | xpmask <xpmask>}]\n"
          "  Check atoms in <mask> for atomic overlaps less than <cut> (default 0.8 Ang)\n"
          "  and unusual bond lengths greater than equilibrium length + <offset>\n"
          "  (default 1.15 Ang) or less than equilibrium length - <minoffset>\n"
          "  (default 0.50 Ang). If 'around' is specified, check between atoms in\n"
          "  <mask1> and <mask2>. If the frame has box info and 'around' is not\n"
          "  specified a pair list will be used to speed up the calculation. The\n"
          "  cutoff used to build the pair list can be adjusted with 'plcut'\n"
          "  (default 4.0 Ang. or <cut>, whichever is greater).\n"
          "  Check for ring-bond intersections unless 'noringcheck' is specified.\n"
          "  Ring-bond intersections are detected when the ring center to bond\n"
          "  center is less than <rsdist>, or less than <ringdcut> and the angle\n"
          "  between the bond vector and ring perpendicular vector is less than\n"
          "  <ringacut>.\n"
          "  By default, extra points will be ignored as determined by 'xpmask'; the\n"
          "  default mask to select extra points is '&@\\XP'. If 'checkxp' is specified\n",
          "  extra points will also be checked.\n"
          "  Warnings will go to the file specified by 'reportfile', STDOUT,\n"
          "  or will be suppressed if 'silent' is specified. If 'skipbadframes'\n"
          "  is specified, subsequent Actions will be skipped if any problems\n"
          "  are detected.\n");
}

// Action_CheckStructure::Init()
Action::RetType Action_CheckStructure::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  // Get Keywords
  std::string around = actionArgs.GetStringKey("around");
  if (!actionArgs.hasKey("silent"))
    outfile_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("reportfile"),
                                         "Structure check", DataFileList::TEXT, true);
  else
    outfile_ = 0;
  double nonbondcut = actionArgs.getKeyDouble("cut",0.8);
  double plcut      = actionArgs.getKeyDouble("plcut", -1);
  //  actionArgs.getKeyDouble("plcut", std::max(8.0, nonbondcut))
  bool checkxp = false;
  std::string xpmask = "";
  if (actionArgs.hasKey("checkxp"))
    checkxp = true;
  else
    xpmask = actionArgs.GetStringKey("xpmask");
  // Structure checker setup
  int err = check_.SetOptions(
    !(actionArgs.hasKey("noimage")),
    !actionArgs.hasKey("nobondcheck"),
    (outfile_ != 0), // saveProblems
    checkxp,
    debugIn,
    actionArgs.GetMaskNext(),
    around,
    xpmask,
    nonbondcut,
    actionArgs.getKeyDouble("offset",1.15),
    actionArgs.getKeyDouble("minoffset", 0.5),
    plcut,
    !actionArgs.hasKey("noringcheck"),
    actionArgs.getKeyDouble("ringshortdist", 0), // 0 means use StructureCheck default
    actionArgs.getKeyDouble("ringdcut", 0),      // 0 means use StructureCheck default
    actionArgs.getKeyDouble("ringacut", 0)       // 0 means use StructureCheck default
  );
  if (err != 0) return Action::ERR;
  // Remaining keywords
  skipBadFrames_ = actionArgs.hasKey("skipbadframes");
  DataFile* dfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  num_problems_ = init.DSL().AddSet( DataSet::INTEGER, actionArgs.GetStringNext(), "CHECK" );
  if (num_problems_ == 0) return Action::ERR;
  if (dfile != 0) dfile->AddDataSet( num_problems_ );
# ifdef MPI
  idx_ = 0;
  if (outfile_ != 0) {
    MetaData md(num_problems_->Meta().Name(), "frame");
    ds_fn_ = init.DSL().AddSet(DataSet::INTEGER, md);
    md.SetAspect("type");
    ds_pt_ = init.DSL().AddSet(DataSet::INTEGER, md);
    md.SetAspect("a1");
    ds_a1_ = init.DSL().AddSet(DataSet::INTEGER, md);
    md.SetAspect("n1");
    ds_n1_ = init.DSL().AddSet(DataSet::STRING, md);
    md.SetAspect("a2");
    ds_a2_ = init.DSL().AddSet(DataSet::INTEGER, md);
    md.SetAspect("n2");
    ds_n2_ = init.DSL().AddSet(DataSet::STRING, md);
    md.SetAspect("dist");
    ds_d_  = init.DSL().AddSet(DataSet::DOUBLE, md);
    if (ds_fn_ == 0 || ds_pt_ == 0 || ds_a1_ == 0 || ds_n1_ == 0 ||
        ds_a2_ == 0 || ds_n2_ == 0 || ds_d_ == 0)
      return Action::ERR;
  }
# endif

  mprintf("    CHECKSTRUCTURE: Checking atoms in mask '%s'", check_.Mask1().MaskString());
  if (check_.Mask2().MaskStringSet())
    mprintf(" around mask '%s'", check_.Mask2().MaskString());
  if (!check_.ImageOpt().UseImage())
    mprintf(", imaging off");
  if (outfile_ != 0)
    mprintf(", warnings output to %s", outfile_->Filename().full());
  else
    mprintf(", warnings suppressed");
  mprintf(".\n");
  mprintf("\tNumber of problems in each frame will be saved to set '%s'\n",
          num_problems_->legend());
  if (dfile != 0)
    mprintf("\tNumber of problems each frame will be written to '%s'\n",
            dfile->DataFilename().full());
  if (!check_.CheckBonds())
    mprintf("\tNot checking bond lengths.\n");
  else
    mprintf("\tChecking for bond lengths greater than eq. length + %.2f Ang and\n"
            "\t  less than eq. length - %.2f Ang.\n", check_.BondOffset(),
            check_.BondMinOffset());
  if (!check_.CheckRings())
    mprintf("\tNot checking for bond-ring intersections.\n");
  else {
    mprintf("\tChecking for bond-ring intersections.\n");
    mprintf("\t  Ring center to bond center distance < %g Ang. is considered intersecting.\n",
            check_.RingShortDist());
    mprintf("\t  Ring center to bond center distance < %g Ang. is considered intersecting\n"
            "\t    if the angle between the bond and perpendicular ring vectors is < %g deg.\n",
            check_.RingDist(), check_.RingAngleCut_Deg());
  }
  mprintf("\tChecking for inter-atomic distances < %.2f Ang.\n",
          nonbondcut);
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
  if ( check_.PairListCut() < 0.0 )
    mprintf("\tPair list cutoff will be determined from system size.\n");
  else
    mprintf("\tCutoff for building pair list is %f Ang.\n", check_.PairListCut());
# ifdef _OPENMP
  mprintf("\tParallelizing calculation with %u threads.\n", check_.Nthreads());
# endif
  return Action::OK;
}

// Action_CheckStructure::Setup()
Action::RetType Action_CheckStructure::Setup(ActionSetup& setup) {
  CurrentParm_ = setup.TopAddress();
  if (check_.Debug() > 0) {
    if (setup.CoordInfo().HasBox())
      setup.CoordInfo().TrajBox().PrintInfo();
  }
  if (check_.Setup( setup.Top(), setup.CoordInfo().TrajBox() )) return Action::ERR;
  check_.Mask1().MaskInfo();
  if (check_.Mask2().MaskStringSet())
    check_.Mask2().MaskInfo();
  if (check_.CheckBonds())
    mprintf("\tChecking %u bonds.\n", check_.Nbonds());
  // Print imaging info for this parm
  if (check_.ImageOpt().ImagingEnabled())
    mprintf("\tImaging on.\n");
  else
    mprintf("\timaging off.\n");

  return Action::OK;
}

/** Consolidate problems from different threads if necessary and write out. */
void Action_CheckStructure::WriteProblems(StructureCheck::FmtType ft, int frameNum, Topology const& top) {
# ifdef MPI
  for (StructureCheck::const_iterator p = check_.begin(); p != check_.end(); ++p) {
    int atom1 = p->A1() + 1;
    int atom2 = p->A2() + 1;
    ds_fn_->Add(idx_, &frameNum);
    ds_pt_->Add(idx_, &ft);
    ds_a1_->Add(idx_, &atom1);
    ds_n1_->Add(idx_, top.TruncResAtomName(p->A1()).c_str());
    ds_a2_->Add(idx_, &atom2);
    ds_n2_->Add(idx_, top.TruncResAtomName(p->A2()).c_str());
    ds_d_->Add(idx_,  p->Dptr());
    idx_++;
  }
# else
  check_.WriteProblemsToFile(outfile_, frameNum, top);
# endif
}

// Action_CheckStructure::DoAction()
Action::RetType Action_CheckStructure::DoAction(int frameNum, ActionFrame& frm) {
# ifdef TIMER
  t_total_.Start();
# endif
  int fnum = frm.TrajoutNum() + 1;
# ifdef TIMER
  t_overlap_.Start();
# endif
  int total_problems = check_.CheckOverlaps(frm.Frm());
# ifdef TIMER
  t_overlap_.Stop();
# endif
  if (outfile_ != 0) WriteProblems(StructureCheck::F_ATOM, fnum, *CurrentParm_);
  if (check_.CheckBonds()) {
#   ifdef TIMER
    t_bonds_.Start();
#   endif
    total_problems += check_.CheckBonds(frm.Frm());
#   ifdef TIMER
    t_bonds_.Stop();
#   endif
    if (outfile_ != 0) WriteProblems(StructureCheck::F_BOND, fnum, *CurrentParm_);
#   ifdef TIMER
    t_rings_.Start();
#   endif
    total_problems += check_.CheckRings(frm.Frm()); // FIXME
#   ifdef TIMER
    t_rings_.Stop();
#   endif
    if (outfile_ != 0) WriteProblems(StructureCheck::F_RING, fnum, *CurrentParm_);
  }
  num_problems_->Add( frameNum, &total_problems );
# ifdef TIMER
  t_total_.Stop();
# endif
  if (total_problems > 0 && skipBadFrames_)
    return Action::SUPPRESS_COORD_OUTPUT;
  return Action::OK;
}

#ifdef MPI
int Action_CheckStructure::SyncAction() {
  if (outfile_ == 0) return 0;
  // Get total number of problems
  std::vector<int> rank_frames( trajComm_.Size() );
  trajComm_.GatherMaster( &idx_, 1, MPI_INT, &(rank_frames[0]) );
  for (int rank = 1; rank < trajComm_.Size(); rank++)
    idx_ += rank_frames[ rank ];
  ds_fn_->Sync( idx_, rank_frames, trajComm_ );
  ds_pt_->Sync( idx_, rank_frames, trajComm_ );
  ds_a1_->Sync( idx_, rank_frames, trajComm_ );
  ds_n1_->Sync( idx_, rank_frames, trajComm_ );
  ds_a2_->Sync( idx_, rank_frames, trajComm_ );
  ds_n2_->Sync( idx_, rank_frames, trajComm_ );
  ds_d_->Sync(  idx_, rank_frames, trajComm_ );
  ds_fn_->SetNeedsSync( false );
  ds_pt_->SetNeedsSync( false );
  ds_a1_->SetNeedsSync( false );
  ds_n1_->SetNeedsSync( false );
  ds_a2_->SetNeedsSync( false );
  ds_n2_->SetNeedsSync( false );
  ds_d_->SetNeedsSync(  false );
  return 0;
}
#endif

void Action_CheckStructure::Print() {
# ifdef TIMER
  t_overlap_.WriteTiming(2, "Check overlaps           :", t_total_.Total());
  t_bonds_.WriteTiming  (2, "Check bond lengths       :", t_total_.Total());
  t_rings_.WriteTiming  (2, "Check ring intersections :", t_total_.Total());
  t_total_.WriteTiming(1, "Total check time :");
# endif
# ifdef MPI
  if (outfile_ != 0) {
    mprintf("    CHECKSTRUCTURE: Writing to %s\n", outfile_->Filename().full());
    for (unsigned int idx = 0; idx != ds_fn_->Size(); idx++) {
      // TODO without casts?
      DataSet_integer const& FN = static_cast<DataSet_integer const&>( *ds_fn_ );
      DataSet_integer const& PT = static_cast<DataSet_integer const&>( *ds_pt_ );
      DataSet_integer const& A1 = static_cast<DataSet_integer const&>( *ds_a1_ );
      DataSet_string  const& N1 = static_cast<DataSet_string  const&>( *ds_n1_ );
      DataSet_integer const& A2 = static_cast<DataSet_integer const&>( *ds_a2_ );
      DataSet_string  const& N2 = static_cast<DataSet_string  const&>( *ds_n2_ );
      DataSet_double  const& DV = static_cast<DataSet_double  const&>( *ds_d_  );
      outfile_->Printf(Fmt_[PT[idx]], FN[idx], A1[idx], N1[idx].c_str(),
                       A2[idx], N2[idx].c_str(), DV[idx]);
    }
  }
# endif
}

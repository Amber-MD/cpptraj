#include "Action_Rmsd.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString
#include "DataSet_Mesh.h"

// CONSTRUCTOR
Action_Rmsd::Action_Rmsd() :
  perres_(false),
  perresout_(0),
  perrescenter_(false),
  perresinvert_(false),
  perresavg_(0),
  masterDSL_(0),
  debug_(0),
  fit_(true),
  rotate_(true),
  useMass_(false),
  rmsd_(0),
  rmatrices_(0)
{ }

void Action_Rmsd::Help() const {
  mprintf("\t[<name>] <mask> [<refmask>] [out filename] [nofit | norotate]\n"
          "\t[mass] [savematrices]\n%s"
          "\t[perres perresout <filename> [perresavg <avgfile>]\n"
          "\t [range <resRange>] [refrange <refRange>]\n"
          "\t [perresmask <additional mask>] [perrescenter] [perresinvert]\n",
          ReferenceAction::Help());
  mprintf("  Calculate coordinate root-mean-squared deviation of atoms in <mask>\n");
}

// Action_Rmsd::Init()
/** Called once before traj processing. Set up reference info. */
Action::RetType Action_Rmsd::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  // Check for keywords
  fit_ = !actionArgs.hasKey("nofit");
  if (fit_)
    rotate_ = !actionArgs.hasKey("norotate");
  useMass_ = actionArgs.hasKey("mass");
  DataFile* outfile = init.DFL().AddDataFile(actionArgs.GetStringKey("out"), actionArgs);
  bool saveMatrices = actionArgs.hasKey("savematrices");
  // Reference keywords
  REF_.InitRef(actionArgs, init.DSL(), fit_, useMass_ );
  // Per-res keywords
  perres_ = actionArgs.hasKey("perres");
  if (perres_) {
    perresout_ = init.DFL().AddDataFile( actionArgs.GetStringKey("perresout") );
    perresinvert_ = actionArgs.hasKey("perresinvert");
    TgtRange_.SetRange( actionArgs.GetStringKey("range") );
    RefRange_.SetRange( actionArgs.GetStringKey("refrange") );
    perresmask_ = actionArgs.GetStringKey("perresmask");
    if (perresmask_.empty()) 
      perresmask_.assign("");
    else {
      // If perresmask does not start with ampersand, insert one.
      if (perresmask_[0] != '&')
        perresmask_ = '&' + perresmask_;
    }
    perrescenter_ = actionArgs.hasKey("perrescenter");
    perresavg_ = init.DFL().AddDataFile( actionArgs.GetStringKey("perresavg") );
  }
  // Get the RMS mask string for target
  std::string tMaskExpr = actionArgs.GetMaskNext();
  tgtMask_.SetMaskString(tMaskExpr);
  // Get the RMS mask string for reference
  std::string rMaskExpr = actionArgs.GetMaskNext();
  if (rMaskExpr.empty())
    rMaskExpr = tMaskExpr;
  REF_.SetRefMask( rMaskExpr );

  // Set up the RMSD data set.
  MetaData md( actionArgs.GetStringNext(), MetaData::M_RMS ); 
  rmsd_ = init.DSL().AddSet(DataSet::DOUBLE, md, "RMSD");
  if (rmsd_==0) return Action::ERR;
  // Add dataset to data file list
  if (outfile != 0) outfile->AddDataSet( rmsd_ );
  // Set up rotation matrix data set if specified
  if (saveMatrices) {
    md.SetAspect("RM");
    if (!fit_) {
      mprinterr("Error: Must be fitting in order to save rotation matrices.\n");
      return Action::ERR;
    }
    rmatrices_ = init.DSL().AddSet(DataSet::MAT3X3, md);
    if (rmatrices_ == 0) return Action::ERR;
  }
# ifdef MPI
  if (REF_.SetTrajComm( init.TrajComm() )) return Action::ERR;
# endif
  mprintf("    RMSD: (%s), reference is %s", tgtMask_.MaskString(),
          REF_.RefModeString().c_str());
  if (!fit_)
    mprintf(", no fitting");
  else {
    mprintf(", with fitting");
    if (!rotate_)
      mprintf(" (no rotation)");
  }
  if (useMass_)
    mprintf(", mass-weighted");
  mprintf(".\n");
  if (rmatrices_ != 0)
    mprintf("\tRotation matrices will be saved to set '%s'\n", rmatrices_->legend());
  // Per-residue RMSD info.
  if (perres_) {
    mprintf("          No-fit RMSD will also be calculated for ");
    if (TgtRange_.Empty()) 
      mprintf("each solute residue");
    else
      mprintf("residues %s",TgtRange_.RangeArg());
    if (!RefRange_.Empty())
      mprintf(" (reference residues %s)",RefRange_.RangeArg());
    mprintf(" using mask [:X%s].\n",perresmask_.c_str());
    if (perresout_ != 0)
      mprintf("          Per-residue output file is %s\n",perresout_->DataFilename().base());
    if (perresavg_ != 0)
      mprintf("          Avg per-residue output file is %s\n",perresavg_->DataFilename().base());
    if (perrescenter_)
      mprintf("          perrescenter: Each residue will be centered prior to RMS calc.\n");
    if (perresinvert_)
      mprintf("          perresinvert: Frames will be written in rows instead of columns.\n");
  }
  if (perres_)
    init.DSL().SetDataSetsPending(true);
  masterDSL_ = init.DslPtr();
  return Action::OK;
}

// Action_Rmsd::perResSetup()
/** Perform setup required for per residue rmsd calculation.
  * Need to set up a target mask, reference mask, and dataset for each
  * residue specified in ResRange.
  * NOTE: Residues in the range arguments from user start at 1, internal
  *       res nums start from 0.
  */
int Action_Rmsd::perResSetup(Topology const& currentParm, Topology const& refParm) {
  Range tgt_range; // Selected target residues
  Range ref_range; // Selected reference residues

  // If no target range previously specified do all solute residues
  if (TgtRange_.Empty()) { 
    tgt_range = currentParm.SoluteResidues();
    tgt_range.ShiftBy(1); // To match user range arg which would start from 1
  } else
    tgt_range = TgtRange_;
  // If the reference range is empty, set it to match the target range
  if (RefRange_.Empty()) 
    ref_range = tgt_range;
  else
    ref_range = RefRange_;
  // Check that the number of reference residues matches number of target residues
  if (tgt_range.Size() != ref_range.Size()) {
    mprintf("Warning: Number of target residues %i does not match\n"
            "Warning:   number of reference residues %i.\n",
            tgt_range.Size(), ref_range.Size());
    return 1;
  }

  // Setup a dataset, target mask, and reference mask for each residue.
  int maxNatom = 0;
  Range::const_iterator ref_it = ref_range.begin();
  MetaData md(rmsd_->Meta().Name(), "res");
  md.SetScalarMode( MetaData::M_RMS );
  for (Range::const_iterator tgt_it = tgt_range.begin();
                             tgt_it != tgt_range.end(); ++tgt_it, ++ref_it)
  {
    int tgtRes = *tgt_it;
    int refRes = *ref_it;
    // Check if either the residue num or the reference residue num out of range.
    if ( tgtRes < 1 || tgtRes > currentParm.Nres()) {
      mprintf("Warning: Specified residue # %i is out of range.\n", tgtRes);
      continue;
    }
    if ( refRes < 1 || refRes > refParm.Nres() ) {
      mprintf("Warning: Specified reference residue # %i is out of range.\n", refRes);
      continue;
    }
    // Check if a perResType has been set for this residue # yet.
    perResArray::iterator PerRes;
    for (PerRes = ResidueRMS_.begin(); PerRes != ResidueRMS_.end(); ++PerRes)
      if ( PerRes->data_->Meta().Idx() == tgtRes ) break;
    // If necessary, create perResType for residue
    if (PerRes == ResidueRMS_.end()) {
      perResType p;
      md.SetIdx( tgtRes );
      md.SetLegend( currentParm.TruncResNameNum(tgtRes-1) );
      p.data_ = (DataSet_1D*)masterDSL_->AddSet(DataSet::DOUBLE, md);
      if (p.data_ == 0) {
        mprinterr("Internal Error: Could not set up per residue data set.\n");
        return 2;
      }
      if (perresout_ != 0) perresout_->AddDataSet( p.data_ );
      // Setup mask strings. Note that masks are based off user residue nums
      p.tgtResMask_.SetMaskString(":" + integerToString(tgtRes) + perresmask_);
      p.refResMask_.SetMaskString(":" + integerToString(refRes) + perresmask_);
      ResidueRMS_.push_back( p );
      PerRes = ResidueRMS_.end() - 1;
    }
    PerRes->isActive_ = false;
    // Setup the reference mask
    if (refParm.SetupIntegerMask(PerRes->refResMask_)) {
      mprintf("Warning: Could not setup reference mask for residue %i\n",refRes);
      continue;
    }
    if (PerRes->refResMask_.None()) {
      mprintf("Warning: No atoms selected for reference residue %i\n",refRes);
      continue;
    }
    // Setup the target mask
    if (currentParm.SetupIntegerMask(PerRes->tgtResMask_)) {
      mprintf("Warning: Could not setup target mask for residue %i\n",tgtRes);
      continue;
    }
    if (PerRes->tgtResMask_.None()) {
      mprintf("Warning: No atoms selected for target residue %i\n",tgtRes);
      continue;
    }
    // Check that # atoms in target and reference masks match
    if (PerRes->tgtResMask_.Nselected() != PerRes->refResMask_.Nselected()) {
      mprintf("Warning: Res %i: # atoms in Tgt (%i) != # atoms in Ref (%i)\n",
              tgtRes, PerRes->tgtResMask_.Nselected(), PerRes->refResMask_.Nselected());
      if (debug_ > 0) {
        mprintf("    Target Atoms:\n");
        for (AtomMask::const_iterator t = PerRes->tgtResMask_.begin();
                                      t != PerRes->tgtResMask_.end(); ++t)
          mprintf("\t%s\n", currentParm.AtomMaskName(*t).c_str());
        mprintf("    Ref Atoms:\n");
        for (AtomMask::const_iterator r = PerRes->refResMask_.begin();
                                      r != PerRes->refResMask_.end(); ++r)
          mprintf("\t%s\n", refParm.AtomMaskName(*r).c_str());
      }
      continue;
    }
    if ( PerRes->tgtResMask_.Nselected() > maxNatom ) maxNatom = PerRes->tgtResMask_.Nselected();
    // Indicate that these masks were properly set up
    PerRes->isActive_ = true;
  }
  mprintf("\tMax # selected atoms in residues: %i\n", maxNatom);

  // Allocate memory for target and reference residue frames.
  // Although initial masses are wrong this is OK since the number of atoms 
  // and masses will be assigned when residue RMSD is actually being calcd.
  if (maxNatom > 0) {
    std::vector<Atom> temp( maxNatom );
    ResTgtFrame_.SetupFrameM( temp );
    ResRefFrame_.SetupFrameM( temp );
  } else {
    mprintf("Warning: No residues selected for per-residue calculation.\n");
    return 1;
  } 
    
  return 0;
}

// Action_Rmsd::Setup()
/** Called every time the trajectory changes. Set up FrameMask for the new 
  * parmtop and allocate space for selected atoms from the Frame.
  */
Action::RetType Action_Rmsd::Setup(ActionSetup& setup) {
  // Target setup
  if ( setup.Top().SetupIntegerMask( tgtMask_ ) ) return Action::ERR;
  mprintf("\tTarget mask:");
  tgtMask_.BriefMaskInfo();
  mprintf("\n");
  if ( tgtMask_.None() ) {
    mprintf("Warning: No atoms in mask '%s'.\n", tgtMask_.MaskString());
    return Action::SKIP;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  tgtFrame_.SetupFrameFromMask(tgtMask_, setup.Top().Atoms());
  // Reference setup
  if (REF_.SetupRef(setup.Top(), tgtMask_.Nselected()))
    return Action::SKIP;
 
  // Per residue rmsd setup
  if (perres_) {
    Topology* RefParm = REF_.RefCrdTopPtr();
    // If RefParm is still NULL probably 'first', set now.
    if (RefParm == 0)
      RefParm = setup.TopAddress();
    int err = perResSetup(setup.Top(), *RefParm);
    if      (err == 1) return Action::SKIP;
    else if (err == 2) return Action::ERR;
  }

  // Warn if PBC and rotating
  if (rotate_ && setup.CoordInfo().TrajBox().Type() != Box::NOBOX) {
    mprintf("Warning: Coordinates are being rotated and box coordinates are present.\n"
            "Warning: Unit cell vectors are NOT rotated; imaging will not be possible\n"
            "Warning:  after the RMS-fit is performed.\n");
  }

  return Action::OK;
}

// Action_Rmsd::DoAction()
Action::RetType Action_Rmsd::DoAction(int frameNum, ActionFrame& frm) {
  // Perform any needed reference actions
  REF_.ActionRef( frm.TrajoutNum(), frm.Frm() );
  // Calculate RMSD
  double rmsdval;
  Action::RetType err;
  // Set selected frame atoms. Masses have already been set.
  tgtFrame_.SetCoordinates(frm.Frm(), tgtMask_);
  if (!fit_) {
    rmsdval = tgtFrame_.RMSD_NoFit(REF_.SelectedRef(), useMass_);
    err = Action::OK;
  } else {
    rmsdval = tgtFrame_.RMSD_CenteredRef(REF_.SelectedRef(), rot_, tgtTrans_, useMass_);
    if (rmatrices_ != 0) rmatrices_->Add(frameNum, rot_.Dptr());
    if (rotate_)
      frm.ModifyFrm().Trans_Rot_Trans(tgtTrans_, rot_, REF_.RefTrans());
    else {
      tgtTrans_ += REF_.RefTrans();
      frm.ModifyFrm().Translate(tgtTrans_);
    }
    err = Action::MODIFY_COORDS;
  }
  rmsd_->Add(frameNum, &rmsdval);

  // ---=== Per Residue RMSD ===---
  // Set reference and selected frame for each residue using the previously
  // set-up masks in refResMask and tgtResMask. Use SetFrame instead
  // of SetCoordinates since each residue can be a different size.
  if (perres_) {
    for (perResArray::const_iterator PerRes = ResidueRMS_.begin();
                                     PerRes != ResidueRMS_.end(); ++PerRes)
    {
      if ( PerRes->isActive_ ) {
        ResRefFrame_.SetFrame(REF_.CurrentReference(), PerRes->refResMask_);
        ResTgtFrame_.SetFrame(frm.Frm(),               PerRes->tgtResMask_);
        if (perrescenter_) {
          ResTgtFrame_.CenterOnOrigin( useMass_ );
          ResRefFrame_.CenterOnOrigin( useMass_ );
        }
        double R = ResTgtFrame_.RMSD_NoFit(ResRefFrame_, useMass_);
        PerRes->data_->Add(frameNum, &R);
      }
    }
  }
  REF_.PreviousRef( frm.Frm() );

  return err;
}

// Action_Rmsd::Print()
/** For per-residue RMSD only. Setup output
  * file options. Calculate averages if requested.
  */
void Action_Rmsd::Print() {
  if (!perres_ || ResidueRMS_.empty()) return;
  // Per-residue output file
  if (perresout_ != 0) {
    // Set output file to be inverted if requested
    if (perresinvert_) 
      perresout_->ProcessArgs("invert");
    mprintf("    RMSD: Per-residue: Writing data for %zu residues to %s\n",
            ResidueRMS_.size(), perresout_->DataFilename().full());
  }

  // Average
  if (perresavg_ != 0) {
    // Use the per residue rmsd dataset list to add one more for averaging
    DataSet_Mesh* PerResAvg = (DataSet_Mesh*)masterDSL_->AddSet(DataSet::XYMESH, 
                                                                MetaData(rmsd_->Meta().Name(),
                                                                         "Avg"));
    PerResAvg->ModifyDim(Dimension::X).SetLabel("Residue");
    // another for stdev
    DataSet_Mesh* PerResStdev = (DataSet_Mesh*)masterDSL_->AddSet(DataSet::XYMESH, 
                                                                  MetaData(rmsd_->Meta().Name(),
                                                                           "Stdev"));
    PerResStdev->ModifyDim(Dimension::X).SetLabel("Residue");
#   ifdef MPI
    PerResAvg->SetNeedsSync( false );
    PerResStdev->SetNeedsSync( false );
#   endif
    // Add the average and stdev datasets to the master datafile list
    perresavg_->AddDataSet(PerResAvg);
    perresavg_->AddDataSet(PerResStdev);
    // For each residue, get the average rmsd
    double stdev = 0;
    double avg = 0;
    for (perResArray::const_iterator PerRes = ResidueRMS_.begin();
                                     PerRes != ResidueRMS_.end(); ++PerRes)
    {
      avg = PerRes->data_->Avg( stdev );
      double pridx = (double)PerRes->data_->Meta().Idx();
      PerResAvg->AddXY(pridx, avg);
      PerResStdev->AddXY(pridx, stdev);
    }
  }
}

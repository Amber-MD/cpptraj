// RMSD
#include "Action_Rmsd.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString
#include "DS_Math.h" // Avg

// TODO: Make all Frames non-pointers

// CONSTRUCTOR
Action_Rmsd::Action_Rmsd() :
  perres_(false),
  NumResidues_(0),
  perresout_(0),
  perrescenter_(false),
  perresinvert_(false),
  perresavg_(0),
  ResFrame_(0),
  ResRefFrame_(0),
  nofit_(false),
  rotate_(true),
  useMass_(false),
  rmsd_(0),
  masterDSL_(0),
  refmode_(UNKNOWN_REF),
  RefParm_(0)
{ }

// DESTRUCTOR
Action_Rmsd::~Action_Rmsd() {
  //mprinterr("RMSD DESTRUCTOR\n");
  if (ResFrame_!=0) delete ResFrame_;
  if (ResRefFrame_!=0) delete ResRefFrame_;
}

// Action_Rmsd::resizeResMasks()
/** For perres rmsd. If the current number of residues is greater than
  * the size of the residue mask lists, allocate as many extra masks
  * as needed. 
  */
void Action_Rmsd::resizeResMasks() {
  AtomMask Blank;
  if (NumResidues_ > (int)tgtResMask_.size()) {
    tgtResMask_.resize(NumResidues_, Blank);
    refResMask_.resize(NumResidues_, Blank);
  }
} 

/** Setup RefMask based on given Topology. Allocate space for selected
  * reference atoms. Save the reference parm as RefParm for use with
  * PerResSetup.
  */
int Action_Rmsd::SetRefMask( Topology* RefParmIn ) {
  // If a reference parm is already set, exit.
  if (RefParm_!=0) return 0;
  if (RefParmIn == 0) return 1;
  RefParm_ = RefParmIn;
  if (RefParm_->SetupIntegerMask( RefMask_ )) return 1;
  if (RefMask_.None()) {
    mprinterr("Error: rmsd: No reference atoms selected for parm %s, [%s]\n",
              RefParm_->c_str(), RefMask_.MaskString());
    return 1;
  }
  SelectedRef_.SetupFrameFromMask( RefMask_, RefParm_->Atoms() );
  return 0;
}

/** Setup selected reference coordinates based on given frame
  * and RefMask. If fitting, the reference is pre-centered
  * and the translation back to the original reference point
  * is stored.
  */
void Action_Rmsd::SetRefStructure( Frame& frameIn ) {
  RefFrame_ = frameIn;
  SelectedRef_.SetCoordinates( RefFrame_, RefMask_ );
  if (!nofit_)
    refTrans_ = SelectedRef_.CenterOnOrigin( useMass_ );
}

void Action_Rmsd::Help() {
  mprintf("rmsd [<name>] <mask> [<refmask>] [out filename] [nofit | norotate] [mass]\n");
  mprintf("     [ first | ref <filename> | refindex <#> |\n");
  mprintf("     reftraj <filename> [parm <parmname> | parmindex <#>] ]\n");
  mprintf("     [ perres perresout <filename> [range <res range>] [refrange <ref res range>]\n");
  mprintf("     [perresmask <addtl mask>] [perresinvert] [perrescenter] perresavg <pravg> ]\n");
}

// Action_Rmsd::init()
/** Called once before traj processing. Set up reference info. */
Action::RetType Action_Rmsd::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  std::string reftrajname;
  Topology* RefParm = 0;
  ReferenceFrame REF;
  // Check for keywords
  nofit_ = actionArgs.hasKey("nofit");
  if (!nofit_)
    rotate_ = !actionArgs.hasKey("norotate");
  useMass_ = actionArgs.hasKey("mass");
  DataFile* outfile = DFL->AddDataFile(actionArgs.GetStringKey("out"), actionArgs);
  // Reference keywords
  refmode_ = UNKNOWN_REF;
  if ( actionArgs.hasKey("first") ) {
    refmode_ = FIRST;
  } else {
    REF = FL->GetFrame( actionArgs );
    if (REF.error()) return Action::ERR;
    if (REF.empty()) {
      reftrajname = actionArgs.GetStringKey("reftraj");
      if (!reftrajname.empty()) {
        RefParm = PFL->GetParm( actionArgs );
        refmode_ = REFTRAJ;
      } else {
        // No reference keywords specified. Default to first.
        mprintf("Warning: rmsd: No reference structure given. Defaulting to first.\n");
        refmode_ = FIRST;
      }
    } else
      refmode_ = REFFRAME;
  }
  // Per-res keywords
  perres_ = actionArgs.hasKey("perres");
  if (perres_) {
    perresout_ = DFL->AddDataFile( actionArgs.GetStringKey("perresout") );
    perresinvert_ = actionArgs.hasKey("perresinvert");
    ResRange_.SetRange( actionArgs.GetStringKey("range") );
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
    perresavg_ = DFL->AddDataFile( actionArgs.GetStringKey("perresavg") );
  }
  // Get the RMS mask string for target 
  std::string mask0 = actionArgs.GetMaskNext();
  FrameMask_.SetMaskString(mask0);
  // Get the RMS mask string for reference
  std::string mask1 = actionArgs.GetMaskNext();
  if (mask1.empty())
    mask1 = mask0;
  RefMask_.SetMaskString( mask1 );

  // Initialize reference if not 'first'.
  if (refmode_ != FIRST) {
    if ( !reftrajname.empty() ) {
      // Reference trajectory
      if (RefParm == 0) {
        mprinterr("Error: rmsd: Could not get parm for reftraj %s\n", reftrajname.c_str());
        return Action::ERR;
      }
      if (SetRefMask( RefParm )!=0) return Action::ERR;
      // Attempt to open reference traj.
      if (RefTraj_.SetupTrajRead( reftrajname, &actionArgs, RefParm)) {
        mprinterr("Error: rmsd: Could not set up reftraj %s\n", reftrajname.c_str());
        return Action::ERR;
      }
      RefFrame_.SetupFrameV(RefParm->Atoms(), RefTraj_.HasVelocity());
      if (RefTraj_.BeginTraj(false)) {
        mprinterr("Error: rmsd: Could not open reftraj %s\n", reftrajname.c_str());
        return Action::ERR;
      }
    } else {
      // Reference Frame
      if (SetRefMask( REF.Parm() ) != 0) return Action::ERR;
      SetRefStructure( *(REF.Coord()) );
    } 
  }

  // Set up the RMSD data set. 
  rmsd_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"RMSD");
  if (rmsd_==0) return Action::ERR;
  rmsd_->SetScalar( DataSet::M_RMS );
  // Add dataset to data file list
  if (outfile != 0) outfile->AddSet( rmsd_ );

  mprintf("    RMSD: (%s), reference is",FrameMask_.MaskString());
  if (refmode_ == FIRST)
    mprintf(" first frame");
  else if (refmode_==REFTRAJ)
    mprintf(" trajectory %s",RefTraj_.FullTrajStr());
  else // REFFRAME
    mprintf (" reference frame %s", REF.FrameName()); 
  mprintf(" (%s)",RefMask_.MaskString());
  if (nofit_)
    mprintf(", no fitting");
  else {
    mprintf(", with fitting");
    if (!rotate_)
      mprintf(" (no rotation)");
  }
  if (useMass_) 
    mprintf(", mass-weighted");
  mprintf(".\n");
  // Per-residue RMSD info.
  if (perres_) {
    mprintf("          No-fit RMSD will also be calculated for ");
    if (ResRange_.Empty()) 
      mprintf("each solute residue");
    else
      mprintf("residues %s",ResRange_.RangeArg());
    if (!RefRange_.Empty())
      mprintf(" (reference residues %s)",RefRange_.RangeArg());
    mprintf(" using mask [:X%s].\n",perresmask_.c_str());
    if (perresout_ != 0)
      mprintf("          Per-residue output file is %s\n",perresout_->Filename());
    if (perresavg_ != 0)
      mprintf("          Avg per-residue output file is %s\n",perresavg_->Filename());
    if (perrescenter_)
      mprintf("          perrescenter: Each residue will be centered prior to RMS calc.\n");
    if (perresinvert_)
      mprintf("          perresinvert: Frames will be written in rows instead of columns.\n");
  }
  masterDSL_ = DSL;
  return Action::OK;
}

// Action_Rmsd::perResSetup()
/** Perform setup required for per residue rmsd calculation.
  * Need to set up a target mask, reference mask, and dataset for each
  * residue specified in ResRange.
  * NOTE: Residues in the range arguments from user start at 1, internal
  *       res nums start from 0.
  */
int Action_Rmsd::perResSetup(Topology* currentParm, Topology* RefParm) {
  Range tgt_range;
  Range ref_range;

  NumResidues_ = currentParm->FinalSoluteRes();

  // If no target range previously specified do all solute residues
  if (ResRange_.Empty()) 
    tgt_range.SetRange(1, NumResidues_+1);
  else
    tgt_range = ResRange_;

  // If the reference range is empty, set it to match the target range
  if (RefRange_.Empty()) 
    ref_range = tgt_range;
  else
    ref_range = RefRange_;

  // Check that the number of reference residues matches number of target residues
  if (tgt_range.Size() != ref_range.Size()) {
    mprintf("Warning: RMSD: perres: Number of residues %i does not match\n",tgt_range.Size());
    mprintf("Warning:       number of reference residues %i.\n",ref_range.Size());
    return 1;
  }

  // Setup a dataset, target mask, and reference mask, for each residue.
  // Since we will only calculate per res rmsd for residues that can be
  // successfully set up, keep track of that as well.
  //mprinterr("DEBUG: Setting up %i masks and data for %s\n",nres,currentParm->parmName);
  resizeResMasks();
  resIsActive_.reserve(NumResidues_);
  resIsActive_.assign(NumResidues_, false);
  int N = -1; // Set to -1 since increment is at top of loop
  Range::const_iterator ref_it = ref_range.begin();
  for (Range::const_iterator tgt_it = tgt_range.begin();
                             tgt_it != tgt_range.end(); ++tgt_it)
  {
    int tgtRes = *tgt_it;
    int refRes = *ref_it;
    ++ref_it;
    // Check if either the residue num or the reference residue num out of range.
    if ( tgtRes < 1 || tgtRes > NumResidues_) {
      mprintf("Warning: Rmsd: perres: Specified residue # %i is out of range.\n",
              tgtRes);
      continue;
    }
    if ( refRes < 1 || refRes > NumResidues_ ) {
      mprintf("Warning: Rmsd: perres: Specified reference residue # %i is out of range.\n",
              refRes);
      continue;
    }
    ++N;
    // Create dataset for res - if already present this returns null 
    DataSet* prDataSet = masterDSL_->AddSetIdxAspect( DataSet::DOUBLE, rmsd_->Name(), tgtRes, "res");
    prDataSet->SetLegend( currentParm->TruncResNameNum(tgtRes-1) );
    PerResRMSD_.push_back( prDataSet );

    // Setup mask strings. Note that masks are based off user residue nums
    std::string tgtArg = ":" + integerToString(tgtRes) + perresmask_;
    tgtResMask_[N].SetMaskString(tgtArg.c_str());
    std::string refArg = ":" + integerToString(refRes) + perresmask_;
    refResMask_[N].SetMaskString(refArg.c_str());
    //mprintf("DEBUG: RMSD: PerRes: Mask %s RefMask %s\n",tgtArg,refArg);

    // Setup the reference mask
    if (RefParm->SetupIntegerMask(refResMask_[N])) {
      mprintf("      perres: Could not setup reference mask for residue %i\n",refRes);
      continue;
    }
    if (refResMask_[N].None()) {
      mprintf("      perres: No atoms selected for reference residue %i\n",refRes);
      continue;
    }

    // Setup the target mask
    if (currentParm->SetupIntegerMask(tgtResMask_[N])) {
      mprintf("      perres: Could not setup target mask for residue %i\n",tgtRes);
      continue;
    }
    if (tgtResMask_[N].None()) {
      mprintf("      perres: No atoms selected for target residue %i\n",tgtRes);
      continue;
    }

    // Check that # atoms in target and reference masks match
    if (tgtResMask_[N].Nselected() != refResMask_[N].Nselected()) {
      mprintf("      perres: Res %i: # atoms in Tgt [%i] != # atoms in Ref [%i]\n",
              tgtRes,tgtResMask_[N].Nselected(),refResMask_[N].Nselected());
      continue;
    }

    // Indicate that these masks were properly set up
    resIsActive_[N]=true;
  }   

  // Allocate memory for residue frame and residue reference frame. The size 
  // of each Frame is initially allocated to the maximum number of atoms.
  // Although initial masses are wrong this is ok since the number of atoms 
  // and masses will change when residue RMSD is actually being calcd.
  if (ResRefFrame_!=0) delete ResRefFrame_;
  ResRefFrame_ = new Frame( RefParm->Atoms() );
  //ResRefFrame->Info("ResRefFrame");
  if (ResFrame_!=0) delete ResFrame_;
  ResFrame_ = new Frame( currentParm->Atoms() );
  //ResFrame->Info("ResFrame");

  return 0;
}

// Action_Rmsd::setup()
/** Called every time the trajectory changes. Set up FrameMask for the new 
  * parmtop and allocate space for selected atoms from the Frame.
  */
Action::RetType Action_Rmsd::Setup(Topology* currentParm, Topology** parmAddress) {
  if ( currentParm->SetupIntegerMask( FrameMask_ ) ) return Action::ERR;
  FrameMask_.MaskInfo();
  if ( FrameMask_.None() ) {
    mprintf("Warning: rmsd: No atoms in mask.\n");
    return Action::ERR;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  SelectedFrame_.SetupFrameFromMask(FrameMask_, currentParm->Atoms());

  // Reference setup if 'first'
  if (refmode_ == FIRST) {
    if ( SetRefMask( currentParm )!=0 ) return Action::ERR;
  }
  
  // Check that num atoms in frame mask from this parm match ref parm mask
  if ( RefMask_.Nselected() != FrameMask_.Nselected() ) {
    mprintf("Warning: Number of atoms in RMS mask (%i) does not equal number of\n",
              FrameMask_.Nselected());
    mprintf("Warning: atoms in reference mask (%i).\n",RefMask_.Nselected());
    return Action::ERR;
  }

  // Per residue rmsd setup
  if (perres_) { 
    if (perResSetup(currentParm, RefParm_)) return Action::ERR;
  }

  return Action::OK;
}

// Action_Rmsd::action()
/** Called every time a frame is read in. Calc RMSD. If not first and not
  * RefTraj SetRefStructure has already been called. When fitting, 
  * SetRefStructure pre-centers the reference coordinates at the origin
  * and puts the translation from origin to reference in Trans[3-5]. 
  */
Action::RetType Action_Rmsd::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  double R;
  Matrix_3x3 U;

  // Perform any needed reference actions
  if (refmode_ == FIRST) {
    SetRefStructure( *currentFrame );
    refmode_ = REFFRAME;
  } else if (refmode_ == REFTRAJ) {
    RefTraj_.GetNextFrame( RefFrame_ );
    SelectedRef_.SetCoordinates(RefFrame_, RefMask_);
    if (!nofit_)
      refTrans_ = SelectedRef_.CenterOnOrigin(useMass_);
  }

  // Set selected frame atoms. Masses have already been set.
  SelectedFrame_.SetCoordinates(*currentFrame, FrameMask_);

  // DEBUG
/*  mprintf("  DEBUG: RMSD: First atom coord in SelectedFrame is : "); 
  SelectedFrame.printAtomCoord(0);
  mprintf("  DEBUG: RMSD: First atom coord in SelectedRef is : ");
  SelectedRef.printAtomCoord(0);
*/

  if (nofit_) {
    R = SelectedFrame_.RMSD_NoFit(SelectedRef_, useMass_);
  } else {
    R = SelectedFrame_.RMSD_CenteredRef(SelectedRef_, U, Trans_, useMass_);
    if (rotate_)
      currentFrame->Trans_Rot_Trans(Trans_, U, refTrans_);
    else {
      Trans_ += refTrans_;
      currentFrame->Translate(Trans_);
    }
  }

  rmsd_->Add(frameNum, &R);

  // ---=== Per Residue RMSD ===---
  // Set reference and selected frame for each residue using the previously
  // set-up masks in refResMask and tgtResMask. Use SetFrame instead
  // of SetCoordinates since each residue can be a different size.
  if (perres_) {
    for (int N=0; N < NumResidues_; ++N) {
      if (!resIsActive_[N]) {
        //mprintf("DEBUG:           [%4i] Not Active.\n",N);
        continue;
      }
      ResRefFrame_->SetFrame(RefFrame_, refResMask_[N]);
      ResFrame_->SetFrame(*currentFrame, tgtResMask_[N]);
      if (perrescenter_) {
        ResFrame_->CenterOnOrigin(false);
        ResRefFrame_->CenterOnOrigin(false);
      }
      R = ResFrame_->RMSD_NoFit(*ResRefFrame_, useMass_);
      //mprintf("DEBUG:           [%4i] Res [%s] nofit RMSD to [%s] = %lf\n",N,
      //        tgtResMask[N]->MaskString(),refResMask[N]->MaskString(),R);
      PerResRMSD_[N]->Add(frameNum, &R);
    }
  }

  return Action::OK;
}

// Action_Rmsd::print()
/** For per-residue RMSD only. Sync the per-residue RMSD data set since
  * it is not part of the master DataSetList in Cpptraj. Setup output
  * file options. Calculate averages if requested.
  */
void Action_Rmsd::Print() {
  if (!perres_ || PerResRMSD_.empty()) return;
  // Per-residue output
  if (perresout_ != 0) {
    // Add data sets to perresout
    for (std::vector<DataSet*>::iterator set = PerResRMSD_.begin();
                                         set != PerResRMSD_.end(); ++set)
      perresout_->AddSet(*set);
    // Set output file to be inverted if requested
    if (perresinvert_) 
      perresout_->ProcessArgs("invert");
    mprintf("    RMSD: Per-residue: Writing data for %zu residues to %s\n",
            PerResRMSD_.size(), perresout_->Filename());
  }

  // Average
  if (perresavg_ != 0) {
    int Nperres = (int)PerResRMSD_.size();
    // Use the per residue rmsd dataset list to add one more for averaging
    DataSet* PerResAvg = masterDSL_->AddSetAspect(DataSet::DOUBLE, rmsd_->Name(), "Avg");
    // another for stdev
    DataSet* PerResStdev = masterDSL_->AddSetAspect(DataSet::DOUBLE, rmsd_->Name(), "Stdev");
    // Add the average and stdev datasets to the master datafile list
    perresavg_->AddSet(PerResAvg);
    perresavg_->AddSet(PerResStdev);
    perresavg_->ProcessArgs("xlabel Residue");
    // For each residue, get the average rmsd
    double stdev = 0;
    double avg = 0;
    for (int pridx = 0; pridx < Nperres; pridx++) {
      avg = DS_Math::Avg(*PerResRMSD_[pridx], &stdev );
      int dsidx = PerResRMSD_[pridx]->Idx() - 1;
      PerResAvg->Add(dsidx, &avg);
      PerResStdev->Add(dsidx,&stdev);
    }
  }
}

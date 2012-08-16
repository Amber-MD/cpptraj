// RMSD
#include "Action_Rmsd.h"
#include "CpptrajStdio.h"

// TODO: Make all Frames non-pointers

// CONSTRUCTOR
Action_Rmsd::Action_Rmsd() :
  perres_(false),
  NumResidues_(0),
  perresout_(NULL),
  perrescenter_(false),
  perresinvert_(false),
  perresavg_(NULL),
  ResFrame_(NULL),
  ResRefFrame_(NULL),
  nofit_(false),
  rmsd_(NULL),
  refmode_(UNKNOWN_REF),
  RefParm_(NULL)
{
  useMass_=false;
}

// DESTRUCTOR
Action_Rmsd::~Action_Rmsd() {
  //mprinterr("RMSD DESTRUCTOR\n");
  if (ResFrame_!=NULL) delete ResFrame_;
  if (ResRefFrame_!=NULL) delete ResRefFrame_;
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
  if (RefParm_!=NULL) return 0;
  if (RefParmIn == NULL) return 1;
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
    SelectedRef_.CenterReference( Trans_+3, useMass_ );
}

// Action_Rmsd::init()
/** Called once before traj processing. Set up reference info.
  * Expected call: 
  * rmsd <name> <mask> [<refmask>] [out filename] [nofit] [mass]
  *      [ first | ref <filename> | refindex <#> | 
  *        reftraj <filename> [parm <parmname> | parmindex <#>] ] 
  *      [ perres perresout <filename> [range <res range>] [refrange <ref res range>] 
  *        [perresmask <addtl mask>] [perresinvert] [perrescenter] perresavg <pravg> ]
  */
int Action_Rmsd::init( ) {
  std::string refname, reftrajname;
  int refindex = -1;
  Topology* refparm = NULL;

  // Check for other keywords
  nofit_ = actionArgs.hasKey("nofit");
  useMass_ = actionArgs.hasKey("mass");
  char *rmsdFile = actionArgs.getKeyString("out",NULL);
  // Reference keywords
  refmode_ = UNKNOWN_REF;
  if ( actionArgs.hasKey("first") ) {
    refmode_ = FIRST;
  } else {
    refindex = actionArgs.getKeyInt("refindex", -1);
    if (actionArgs.hasKey("reference")) refindex = 0;
    refname = actionArgs.GetStringKey("ref");
    if (refindex==-1 && refname.empty()) {
      reftrajname = actionArgs.GetStringKey("reftraj");
      if (!reftrajname.empty()) {
        refparm = PFL->GetParm( actionArgs );
        refmode_ = REFTRAJ;
      }
    } else
      refmode_ = REF;
  }
  // Per-res keywords
  perres_ = actionArgs.hasKey("perres");
  if (perres_) {
    perresout_ = actionArgs.getKeyString("perresout",NULL);
    perresinvert_ = actionArgs.hasKey("perresinvert");
    ResRange_.SetRange( actionArgs.getKeyString("range",NULL) );
    RefRange_.SetRange( actionArgs.getKeyString("refrange",NULL) );
    perresmask_ = actionArgs.GetStringKey("perresmask");
    if (perresmask_.empty()) perresmask_.assign("");
    perrescenter_ = actionArgs.hasKey("perrescenter");
    perresavg_ = actionArgs.getKeyString("perresavg",NULL);
  }
  // Get the RMS mask string for target 
  char* mask0 = actionArgs.getNextMask();
  FrameMask_.SetMaskString(mask0);
  // Get the RMS mask string for reference
  char* mask1 = actionArgs.getNextMask();
  if (mask1==NULL)
    mask1 = mask0;
  RefMask_.SetMaskString( mask1 );

  // Check that a reference keyword was specified. If not, default to first.
  if (refmode_==UNKNOWN_REF) {
    mprintf("Warning: rmsd: No reference structure given. Defaulting to first.\n");
    refmode_ = FIRST;
  }
  // Initialize reference if not 'first'.
  if (refmode_ != FIRST) {
    if ( !reftrajname.empty() ) {
      // Reference trajectory
      if (refparm == NULL) {
        mprinterr("Error: rmsd: Could not get parm for reftraj %s\n", reftrajname.c_str());
        return 1;
      }
      if (SetRefMask( refparm )!=0) return 1;
      // Attempt to open reference traj.
      if (RefTraj_.SetupRead( reftrajname.c_str(), NULL, refparm)) {
        mprinterr("Error: rmsd: Could not set up reftraj %s\n", reftrajname.c_str());
        return 1;
      }
      RefFrame_.SetupFrameV(refparm->Atoms(), RefTraj_.HasVelocity());
      if (RefTraj_.BeginTraj(false)) {
        mprinterr("Error: rmsd: Could not open reftraj %s\n", reftrajname.c_str());
        return 1;
      }
    } else {
      // Reference by name/tag
      if (!refname.empty())
        refindex = FL->FindName( refname );
      // Get reference by index
      Frame* TempFrame = FL->GetFrame( refindex );
      // Get parm for reference
      refparm = FL->GetFrameParm( refindex );
      if (refparm == NULL) {
        mprinterr("Error: rmsd: Could not get parm for frame %s\n", FL->FrameName(refindex));
        return 1;
      }
      if (SetRefMask( refparm )!=0) return 1;
      SetRefStructure( *TempFrame );
    } 
  }

  // Set up the RMSD data set. 
  rmsd_ = DSL->Add(DataSet::DOUBLE, actionArgs.getNextString(),"RMSD");
  if (rmsd_==NULL) return 1;
  rmsd_->SetScalar( DataSet::M_RMS );
  // Add dataset to data file list
  DFL->Add(rmsdFile, rmsd_);

  //rmsd->Info();
  mprintf("    RMSD: (%s), reference is",FrameMask_.MaskString());
  if (refmode_ == FIRST)
    mprintf(" first frame");
  else if (refmode_==REFTRAJ)
    mprintf(" trajectory %s",RefTraj_.FullTrajStr());
  else {
    mprintf(" reference frame");
    if (!refname.empty())
      mprintf(" %s",refname.c_str());
    else
      mprintf(" index %i", refindex);
  }
  mprintf(" (%s)",RefMask_.MaskString());

  if (nofit_)
    mprintf(", no fitting");
  else
    mprintf(", with fitting");
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
    /*if (perresout_==NULL && perresavg_==NULL) {
      mprinterr("Error: perres specified but no output filename given (perresout | perresavg).\n");
      perres_=false;
      return 1;
    }*/
    if (perresout_!=NULL)
      mprintf("          Per-residue output file is %s\n",perresout_);
    if (perresavg_!=NULL)
      mprintf("          Avg per-residue output file is %s\n",perresavg_);
    if (perrescenter_)
      mprintf("          perrescenter: Each residue will be centered prior to RMS calc.\n");
    if (perresinvert_)
      mprintf("          perresinvert: Frames will be written in rows instead of columns.\n");
  }

  return 0;
}

// Action_Rmsd::perResSetup()
/** Perform setup required for per residue rmsd calculation.
  * Need to set up a target mask, reference mask, and dataset for each
  * residue specified in ResRange.
  * NOTE: Residues in the range arguments from user start at 1, internal
  *       res nums start from 0.
  */
int Action_Rmsd::perResSetup(Topology *RefParm) {
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
    // Create dataset for res - if already present this returns NULL
    DataSet* prDataSet = DSL->AddSetIdxAspect( DataSet::DOUBLE, rmsd_->Name(), tgtRes, "res");
    prDataSet->SetLegend( currentParm->ResNameNum(tgtRes-1) );
    PerResRMSD_.push_back( prDataSet );
    //if (prDataSet != NULL) DFL->Add(perresout_, prDataSet);

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

  // Check pointer to the output file
  /*if (perresout_!=NULL) {
    if (DFL->GetDataFile(perresout_)==NULL) {
      mprinterr("Error: RMSD: Perres output file could not be set up.\n");
      return 1;
    }
  }*/

  // Allocate memory for residue frame and residue reference frame. The size 
  // of each Frame is initially allocated to the maximum number of atoms.
  // Although initial masses are wrong this is ok since the number of atoms 
  // and masses will change when residue RMSD is actually being calcd.
  if (ResRefFrame_!=NULL) delete ResRefFrame_;
  ResRefFrame_ = new Frame( RefParm->Atoms() );
  //ResRefFrame->Info("ResRefFrame");
  if (ResFrame_!=NULL) delete ResFrame_;
  ResFrame_ = new Frame( currentParm->Atoms() );
  //ResFrame->Info("ResFrame");

  return 0;
}

// Action_Rmsd::setup()
/** Called every time the trajectory changes. Set up FrameMask for the new 
  * parmtop and allocate space for selected atoms from the Frame.
  */
int Action_Rmsd::setup() {
  if ( currentParm->SetupIntegerMask( FrameMask_ ) ) return 1;
  if ( FrameMask_.None() ) {
    mprintf("Warning: rmsd: No atoms in mask.\n");
    return 1;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  SelectedFrame_.SetupFrameFromMask(FrameMask_, currentParm->Atoms());

  // Reference setup if 'first'
  if (refmode_ == FIRST) {
    if ( SetRefMask( currentParm )!=0 ) return 1;
  }
  
  // Check that num atoms in frame mask from this parm match ref parm mask
  if ( RefMask_.Nselected() != FrameMask_.Nselected() ) {
    mprintf("Warning: Number of atoms in RMS mask (%i) does not equal number of\n",
              FrameMask_.Nselected());
    mprintf("Warning: atoms in reference mask (%i).\n",RefMask_.Nselected());
    return 1;
  }

  // Per residue rmsd setup
  if (perres_) { 
    if (perResSetup(RefParm_)) return 1;
  }

  mprintf("\t%i atoms selected.\n",FrameMask_.Nselected());

  return 0;
}

// Action_Rmsd::action()
/** Called every time a frame is read in. Calc RMSD. If not first and not
  * RefTraj SetRefStructure has already been called. When fitting, 
  * SetRefStructure pre-centers the reference coordinates at the origin
  * and puts the translation from origin to reference in Trans[3-5]. 
  */
int Action_Rmsd::action() {
  double R, U[9];

  // Perform any needed reference actions
  if (refmode_ == FIRST) {
    SetRefStructure( *currentFrame );
    refmode_ = REF;
  } else if (refmode_ == REFTRAJ) {
    RefTraj_.GetNextFrame( RefFrame_ );
    SelectedRef_.SetCoordinates(RefFrame_, RefMask_);
    if (!nofit_)
      SelectedRef_.CenterReference(Trans_+3, useMass_);
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
    R = SelectedFrame_.RMSD(SelectedRef_, useMass_);
  } else {
    R = SelectedFrame_.RMSD_CenteredRef(SelectedRef_, U, Trans_, useMass_);
    currentFrame->Trans_Rot_Trans(Trans_, U);
  }

  rmsd_->Add(frameNum, &R);

  // ---=== Per Residue RMSD ===---
  // Set reference and selected frame for each residue using the previously
  // set-up masks in refResMask and tgtResMask. Use SetFrameFromMask instead
  // of SetFrameCoordsFromMask since each residue can be a different size.
  if (perres_) {
    for (int N=0; N < NumResidues_; ++N) {
      if (!resIsActive_[N]) {
        //mprintf("DEBUG:           [%4i] Not Active.\n",N);
        continue;
      }
      ResRefFrame_->SetFrame(RefFrame_, refResMask_[N]);
      ResFrame_->SetFrame(*currentFrame, tgtResMask_[N]);
      if (perrescenter_) {
        ResFrame_->ShiftToGeometricCenter( );
        ResRefFrame_->ShiftToGeometricCenter( );
      }
      R = ResFrame_->RMSD(*ResRefFrame_, useMass_);
      //mprintf("DEBUG:           [%4i] Res [%s] nofit RMSD to [%s] = %lf\n",N,
      //        tgtResMask[N]->MaskString(),refResMask[N]->MaskString(),R);
      PerResRMSD_[N]->Add(frameNum, &R);
    }
  }

  return 0;
}

// Action_Rmsd::print()
/** For per-residue RMSD only. Sync the per-residue RMSD data set since
  * it is not part of the master DataSetList in Cpptraj. Setup output
  * file options. Calculate averages if requested.
  */
void Action_Rmsd::print() {
  DataFile *outFile;

  if (!perres_ || PerResRMSD_.empty()) return;
  // Per-residue output
  if (perresout_ != NULL) {
    // Add data sets to perresout
    for (std::vector<DataSet*>::iterator set = PerResRMSD_.begin();
                                         set != PerResRMSD_.end(); ++set)
    {
      outFile = DFL->Add(perresout_, *set);
      if (outFile == NULL) 
        mprinterr("Error adding set %s to file %s\n", (*set)->c_str(), perresout_);
    }
    if (outFile!=NULL) {
      // Set output file to be inverted if requested
      if (perresinvert_) 
        outFile->ProcessArgs("invert");
      mprintf("    RMSD: Per-residue: Writing data for %zu residues to %s\n",
              PerResRMSD_.size(), outFile->Filename());
    }
  }

  // Average
  if (perresavg_ != NULL) {
    int Nperres = (int)PerResRMSD_.size();
    // Use the per residue rmsd dataset list to add one more for averaging
    DataSet *PerResAvg = DSL->AddSetAspect(DataSet::DOUBLE, rmsd_->Name(), "Avg");
    // another for stdev
    DataSet *PerResStdev = DSL->AddSetAspect(DataSet::DOUBLE, rmsd_->Name(), "Stdev");
    // Add the average and stdev datasets to the master datafile list
    outFile = DFL->Add(perresavg_, PerResAvg);
    outFile = DFL->Add(perresavg_, PerResStdev);
    outFile->ProcessArgs("xlabel Residue");
    // For each residue, get the average rmsd
    double stdev = 0;
    double avg = 0;
    for (int pridx = 0; pridx < Nperres; pridx++) {
      avg = PerResRMSD_[pridx]->Avg( &stdev );
      int dsidx = PerResRMSD_[pridx]->Idx() - 1;
      PerResAvg->Add(dsidx, &avg);
      PerResStdev->Add(dsidx,&stdev);
    }
  }
}
 

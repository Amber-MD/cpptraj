// RMSD
//#include <cstdio> // for sprintf
#include "Action_Rmsd.h"
#include "CpptrajStdio.h"
#include "DataSet_double.h" // for SeparateInit

// TODO: Make all Frames non-pointers

// CONSTRUCTOR
Action_Rmsd::Action_Rmsd() :
  perres_(false),
  NumResidues_(0),
  PerResRMSD_(NULL),
  perresout_(NULL),
  perrescenter_(false),
  perresinvert_(false),
  perresavg_(NULL),
  ResFrame_(NULL),
  ResRefFrame_(NULL),
  nofit_(false),
  rmsd_(NULL)
{
  useMass_=false;
}

// DESTRUCTOR
Action_Rmsd::~Action_Rmsd() {
  //mprinterr("RMSD DESTRUCTOR\n");
  if (ResFrame_!=NULL) delete ResFrame_;
  if (ResRefFrame_!=NULL) delete ResRefFrame_;
  if (PerResRMSD_!=NULL) delete PerResRMSD_;
  // If separate, clean up the dataset
  if (isSeparate_) delete rmsd_;
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
  // Check for other keywords
  nofit_ = actionArgs.hasKey("nofit");
  useMass_ = actionArgs.hasKey("mass");
  char *rmsdFile = actionArgs.getKeyString("out",NULL);
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
  char *mask0 = actionArgs.getNextMask();
  FrameMask_.SetMaskString(mask0);

  // Initialize reference. If no reference mask is given mask0 will be used.
  if (RefInit(nofit_, useMass_, mask0, actionArgs, FL, PFL, Trans_+3))
    return 1;

  // Set up the RMSD data set. 
  rmsd_ = DSL->Add(DataSet::DOUBLE, actionArgs.getNextString(),"RMSD");
  if (rmsd_==NULL) return 1;
  rmsd_->SetScalar( DataSet::M_RMS );
  // Add dataset to data file list
  DFL->Add(rmsdFile, rmsd_);

  //rmsd->Info();
  mprintf("    RMSD: (%s), reference is",FrameMask_.MaskString());
  RefInfo();
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
    if (perresout_==NULL && perresavg_==NULL) {
      mprinterr("Error: perres specified but no output filename given (perresout | perresavg).\n");
      perres_=false;
      return 1;
    }
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

// Action_Rmsd::SeparateInit()
/** This routine allows the RMSD action to be initialized outside the main
  * action list so it can be used e.g. by other actions etc. The dataset
  * is allocated locally.
  */
/*int Action_Rmsd::SeparateInit(char *mask0, bool massIn, int debugIn) {
  isSeparate_ = true;
  debug = debugIn;
  useMass_ = massIn;
  // Also set useMassOriginalValue since this is NOT called from 
  // Init.
  useMassOriginalValue_ = useMass_;
  // Only first for reference for now
  SetFirst(nofit, mask0, useMass_); 

  // Set the RMS mask string for target and reference
  FrameMask.SetMaskString(mask0);

  // Set up the RMSD data set. In case the action is being re-initialized,
  // only do this if rmsd is NULL.
  if (rmsd==NULL) {
    rmsd = new DataSet_double();
    if (rmsd->Setup((char*)"RMSD",-1)) return 1;
  }
  return 0;
}*/

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
  int tgtRes, refRes;

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
  if (PerResRMSD_==NULL) PerResRMSD_ = new DataSetList();
  resIsActive_.reserve(NumResidues_);
  resIsActive_.assign(NumResidues_, false);
  int N = -1; // Set to -1 since increment is at top of loop
  tgt_range.Begin();
  ref_range.Begin();
  while (tgt_range.NextInRange(&tgtRes)) {
    ref_range.NextInRange(&refRes);
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
    DataSet* prDataSet = PerResRMSD_->AddMultiN(DataSet::DOUBLE, "", 
                                                currentParm->ResidueName(tgtRes-1),
                                                tgtRes);
    if (prDataSet != NULL) DFL->Add(perresout_, prDataSet);

    // Setup mask strings. Note that masks are based off user residue nums
    std::string tgtArg = ":" + integerToString(tgtRes) + perresmask_;
    //sprintf(tgtArg,":%i%s",tgtRes,perresmask);
    tgtResMask_[N].SetMaskString(tgtArg.c_str());
    std::string refArg = ":" + integerToString(refRes) + perresmask_;
    //sprintf(refArg,":%i%s",refRes,perresmask);
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
  if (perresout_!=NULL) {
    if (DFL->GetDataFile(perresout_)==NULL) {
      mprinterr("Error: RMSD: Perres output file could not be set up.\n");
      return 1;
    }
  }

  // Allocate memory for residue frame and residue reference frame. The size 
  // of each Frame is initially allocated to the maximum number of atoms.
  // Although initial masses are wrong this is ok since the number of atoms 
  // and masses will change when residue RMSD is actually being calcd.
  if (ResRefFrame_!=NULL) delete ResRefFrame_;
  ResRefFrame_ = new Frame( RefParm->FindResidueMaxNatom(), RefParm->Mass() );
  //ResRefFrame->Info("ResRefFrame");
  if (ResFrame_!=NULL) delete ResFrame_;
  ResFrame_ = new Frame( currentParm->FindResidueMaxNatom(), currentParm->Mass() );
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
  SelectedFrame_.SetupFrameFromMask(FrameMask_, currentParm->Mass());

  // Reference setup
  if (RefSetup( currentParm )) return 1;
  
  // Check that num atoms in frame mask from this parm match ref parm mask
  if ( RefNselected() != FrameMask_.Nselected() ) {
    mprintf("Warning: Number of atoms in RMS mask (%i) does not equal number of\n",
              FrameMask_.Nselected());
    mprintf("Warning: atoms in reference mask (%i).\n",RefNselected());
    return 1;
  }

  // Per residue rmsd setup
  if (perres_) { 
    if (this->perResSetup(GetRefParm())) return 1;
  }

  if (!isSeparate_)
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
  RefAction(currentFrame, Trans_+3);

  // Set selected frame atoms. Masses have already been set.
  SelectedFrame_.SetCoordinates(*currentFrame, FrameMask_);

  // DEBUG
/*  mprintf("  DEBUG: RMSD: First atom coord in SelectedFrame is : "); 
  SelectedFrame.printAtomCoord(0);
  mprintf("  DEBUG: RMSD: First atom coord in SelectedRef is : ");
  SelectedRef.printAtomCoord(0);
*/

  if (nofit_) {
    R = SelectedFrame_.RMSD(&SelectedRef_, useMass_);
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
      R = ResFrame_->RMSD(ResRefFrame_, useMass_);
      //mprintf("DEBUG:           [%4i] Res [%s] nofit RMSD to [%s] = %lf\n",N,
      //        tgtResMask[N]->MaskString(),refResMask[N]->MaskString(),R);
      // NOTE: Should check for error on AddData?
      PerResRMSD_->AddData(frameNum, &R, N);
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

  if (!perres_ || PerResRMSD_==NULL) return;
  // Sync dataset list here since it is not part of master dataset list
  PerResRMSD_->Sync();
  // Per-residue output
  outFile = DFL->GetDataFile(perresout_);
  if (outFile!=NULL) {
    // Set output file to be inverted if requested
    if (perresinvert_) 
      outFile->ProcessArgs("invert");
    mprintf("    RMSD: Per-residue: Writing data for %i residues to %s\n",
            PerResRMSD_->Size(), outFile->Filename());
  }

  // Average
  if (perresavg_==NULL) return;
  int Nperres = PerResRMSD_->Size();
  // Use the per residue rmsd dataset list to add one more for averaging
  DataSet *PerResAvg = PerResRMSD_->Add(DataSet::DOUBLE, (char*)"AvgRMSD", "AvgRMSD");
  // another for stdev
  DataSet *PerResStdev = PerResRMSD_->Add(DataSet::DOUBLE, (char*)"Stdev", "Stdev");
  // Add the average and stdev datasets to the master datafile list
  outFile = DFL->Add(perresavg_, PerResAvg);
  outFile = DFL->Add(perresavg_, PerResStdev);
  outFile->ProcessArgs("xlabel Residue");
  // For each residue, get the average rmsd
  double stdev = 0;
  double avg = 0;
  for (int pridx = 0; pridx < Nperres; pridx++) {
    DataSet *tempDS = PerResRMSD_->GetDataSetN(pridx);
    avg = tempDS->Avg(&stdev);
    int dsidx = tempDS->Idx() - 1; // When set up actual resnum is used - change?
    PerResAvg->Add(dsidx, &avg);
    PerResStdev->Add(dsidx,&stdev);
  }
}
 

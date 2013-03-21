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
  RefParm_(0),
  rmsd_(0),
  masterDSL_(0)
{ }

// DESTRUCTOR
Action_Rmsd::~Action_Rmsd() {
  //mprinterr("RMSD DESTRUCTOR\n");
  if (ResFrame_!=0) delete ResFrame_;
  if (ResRefFrame_!=0) delete ResRefFrame_;
}

void Action_Rmsd::Help() {
  mprintf("\t[<name>] <mask> [<refmask>] [out filename] [nofit | norotate] [mass]\n");
  mprintf("\t[ first | ref <filename> | refindex <#> |\n");
  mprintf("\t  reftraj <filename> [parm <parmname> | parmindex <#>] ]\n");
  mprintf("\t[perres perresout <filename> [perresavg <avgfile>]\n");
  mprintf("\t [range <resRange>] [refrange <refRange>]\n");
  mprintf("\t [perresmask <additional mask>] [perrescenter] [perresinvert]\n");
  mprintf("\tCalculate coordinate root-mean-squared deviation of atoms in <mask>\n");
}

// Action_Rmsd::init()
/** Called once before traj processing. Set up reference info. */
Action::RetType Action_Rmsd::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Check for keywords
  GetRmsKeywords( actionArgs );
  DataFile* outfile = DFL->AddDataFile(actionArgs.GetStringKey("out"), actionArgs);
  // Reference keywords
  bool previous = actionArgs.hasKey("previous");
  bool first = actionArgs.hasKey("first");
  ReferenceFrame  REF = FL->GetFrameFromArgs( actionArgs );
  std::string reftrajname = actionArgs.GetStringKey("reftraj");
  RefParm_ = PFL->GetParm( actionArgs );
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
  std::string mask1 = GetRmsMasks(actionArgs); 
  // Initialize reference
  if (InitRef(previous, first, UseMass(), Fit(), reftrajname, REF, RefParm_, mask1,
              actionArgs, "rmsd"))
    return Action::ERR;
  // Set RefParm for perres if not empty
  if (perres_ && RefParm_ == 0 && !REF.empty())
    RefParm_ = REF.Parm();

  // Set up the RMSD data set. 
  rmsd_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"RMSD");
  if (rmsd_==0) return Action::ERR;
  rmsd_->SetScalar( DataSet::M_RMS );
  // Add dataset to data file list
  if (outfile != 0) outfile->AddSet( rmsd_ );

  mprintf("    RMSD: (%s), reference is %s",TgtMask().MaskString(),
          RefModeString());
  PrintRmsStatus();
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
      mprintf("          Per-residue output file is %s\n",perresout_->DataFilename().base());
    if (perresavg_ != 0)
      mprintf("          Avg per-residue output file is %s\n",perresavg_->DataFilename().base());
    if (perrescenter_)
      mprintf("          perrescenter: Each residue will be centered prior to RMS calc.\n");
    if (perresinvert_)
      mprintf("          perresinvert: Frames will be written in rows instead of columns.\n");
  }
  masterDSL_ = DSL;
  return Action::OK;
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
  // Target setup
  if (SetupRmsMask(*currentParm, "rmsd")) return Action::ERR;
  // Reference setup
  if (SetupRef(*currentParm, TgtMask().Nselected(), "rmsd"))
    return Action::ERR;
 
  // Per residue rmsd setup
  if (perres_) {
    // If RefParm is still NULL probably 'first', set now.
    if (RefParm_ == 0)
      RefParm_ = currentParm;
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
  // Perform any needed reference actions
  ActionRef( *currentFrame, Fit(), UseMass() );
  // Calculate RMSD
  double rmsdval = CalcRmsd( *currentFrame, SelectedRef(), RefTrans() );
  rmsd_->Add(frameNum, &rmsdval);

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
      ResRefFrame_->SetFrame(RefFrame(), refResMask_[N]);
      ResFrame_->SetFrame(*currentFrame, tgtResMask_[N]);
      if (perrescenter_) {
        ResFrame_->CenterOnOrigin(false);
        ResRefFrame_->CenterOnOrigin(false);
      }
      double R = ResFrame_->RMSD_NoFit(*ResRefFrame_, UseMass());
      //mprintf("DEBUG:           [%4i] Res [%s] nofit RMSD to [%s] = %lf\n",N,
      //        tgtResMask[N]->MaskString(),refResMask[N]->MaskString(),R);
      PerResRMSD_[N]->Add(frameNum, &R);
    }
  }

  if (Previous())
    SetRefStructure( *currentFrame, Fit(), UseMass() );

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
            PerResRMSD_.size(), perresout_->DataFilename().base());
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

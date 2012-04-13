// RMSD
#include <cstdio> // for sprintf
#include "Action_Rmsd.h"
#include "CpptrajStdio.h"
#include "DataSet_double.h" // for SeparateInit

// TODO: Make all Frames non-pointers

// CONSTRUCTOR
Rmsd::Rmsd() {
  rmsd=NULL;
  PerResRMSD=NULL;
  NumResidues=0;
  perresout=NULL;
  nofit=false;
  useMass_=false;
  perres=false;
  ResFrame=NULL;
  ResRefFrame=NULL;
  perresmask=NULL;
  perrescenter=false;
  perresinvert=false;
  perresavg=NULL;
}

// DESTRUCTOR
Rmsd::~Rmsd() {
  //mprinterr("RMSD DESTRUCTOR\n");
  if (ResFrame!=NULL) delete ResFrame;
  if (ResRefFrame!=NULL) delete ResRefFrame;
  if (PerResRMSD!=NULL) delete PerResRMSD;
  // If separate, clean up the dataset
  if (isSeparate_) delete rmsd;
}

// Rmsd::resizeResMasks()
/** For perres rmsd. If the current number of residues is greater than
  * the size of the residue mask lists, allocate as many extra masks
  * as needed. 
  */
void Rmsd::resizeResMasks() {
  AtomMask Blank;
  if (NumResidues > (int)tgtResMask.size()) {
    tgtResMask.resize(NumResidues, Blank);
    refResMask.resize(NumResidues, Blank);
  }
} 
    
// Rmsd::init()
/** Called once before traj processing. Set up reference info.
  * Expected call: 
  * rmsd <name> <mask> [<refmask>] [out filename] [nofit] [mass]
  *      [ first | ref <filename> | refindex <#> | 
  *        reftraj <filename> [parm <parmname> | parmindex <#>] ] 
  *      [ perres perresout <filename> [range <res range>] [refrange <ref res range>] 
  *        [perresmask <addtl mask>] [perresinvert] [perrescenter] perresavg <pravg> ]
  */
int Rmsd::init( ) {

  // Check for other keywords
  nofit = actionArgs.hasKey("nofit");
  useMass_ = actionArgs.hasKey("mass");
  char *rmsdFile = actionArgs.getKeyString("out",NULL);
  // Per-res keywords
  perres = actionArgs.hasKey("perres");
  if (perres) {
    perresout = actionArgs.getKeyString("perresout",NULL);
    perresinvert = actionArgs.hasKey("perresinvert");
    ResRange.SetRange( actionArgs.getKeyString("range",NULL) );
    RefRange.SetRange( actionArgs.getKeyString("refrange",NULL) );
    perresmask = actionArgs.getKeyString("perresmask",(char*)"");
    perrescenter = actionArgs.hasKey("perrescenter");
    perresavg = actionArgs.getKeyString("perresavg",NULL);
  }
  // Get the RMS mask string for target 
  char *mask0 = actionArgs.getNextMask();
  FrameMask.SetMaskString(mask0);

  // Initialize reference. If no reference mask is given mask0 will be used.
  if (RefInit(nofit, useMass_, mask0, actionArgs, FL, PFL, Trans+3))
    return 1;

  // Set up the RMSD data set. 
  rmsd = DSL->Add(DataSet::DOUBLE, actionArgs.getNextString(),"RMSD");
  if (rmsd==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(rmsdFile,rmsd);

  //rmsd->Info();
  mprintf("    RMSD: (%s), reference is",FrameMask.MaskString());
  RefInfo();
  if (nofit)
    mprintf(", no fitting");
  else
    mprintf(", with fitting");
  if (useMass_) 
    mprintf(", mass-weighted");
  mprintf(".\n");
  // Per-residue RMSD info.
  if (perres) {
    mprintf("          No-fit RMSD will also be calculated for ");
    if (ResRange.Empty()) 
      mprintf("each solute residue");
    else
      mprintf("residues %s",ResRange.RangeArg());
    if (!RefRange.Empty())
      mprintf(" (reference residues %s)",RefRange.RangeArg());
    mprintf(" using mask [:X%s].\n",perresmask);
    if (perresout==NULL && perresavg==NULL) {
      mprintf("Error: perres specified but no output filename given (perresout | perresavg).\n");
      perres=false;
      return 1;
    }
    if (perresout!=NULL)
      mprintf("          Per-residue output file is %s\n",perresout);
    if (perresavg!=NULL)
      mprintf("          Avg per-residue output file is %s\n",perresavg);
    if (perrescenter)
      mprintf("          perrescenter: Each residue will be centered prior to RMS calc.\n");
    if (perresinvert)
      mprintf("          perresinvert: Frames will be written in rows instead of columns.\n");
  }

  return 0;
}

// Rmsd::SeparateInit()
/** This routine allows the RMSD action to be initialized outside the main
  * action list so it can be used e.g. by other actions etc. The dataset
  * is allocated locally.
  */
int Rmsd::SeparateInit(char *mask0, bool massIn, int debugIn) {
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
}

// Rmsd::perResSetup()
/** Perform setup required for per residue rmsd calculation.
  * Need to set up a target mask, reference mask, and dataset for each
  * residue specified in ResRange.
  * NOTE: Residues in the range arguments from user start at 1, internal
  *       res nums start from 0.
  */
int Rmsd::perResSetup(Topology *RefParm) {
  char tgtArg[1024];
  char refArg[1024];
  Range tgt_range;
  Range ref_range;
  int N, tgtRes, refRes;

  NumResidues = currentParm->FinalSoluteRes();

  // If no target range previously specified do all solute residues
  if (ResRange.Empty()) 
    tgt_range.SetRange(1,NumResidues+1);
  else
    tgt_range.SetRange(&ResRange);

  // If the reference range is empty, set it to match the target range
  if (RefRange.Empty()) 
    ref_range.SetRange(&tgt_range);
  else
    ref_range.SetRange(&RefRange);

  // Check that the number of reference residues matches number of target residues
  if (tgt_range.Size() != ref_range.Size()) {
    mprinterr("Error: RMSD: PerRes: Number of residues %i does not match\n",tgt_range.Size());
    mprinterr("       number of reference residues %i.\n",ref_range.Size());
    return 1;
  }

  // Setup a dataset, target mask, and reference mask, for each residue.
  // Since we will only calculate per res rmsd for residues that can be
  // successfully set up, keep track of that as well.
  //mprinterr("DEBUG: Setting up %i masks and data for %s\n",nres,currentParm->parmName);
  resizeResMasks();
  if (PerResRMSD==NULL) PerResRMSD=new DataSetList();
  resIsActive.reserve(NumResidues);
  resIsActive.assign(NumResidues,false);
  N = -1; // Set to -1 since increment is at top of loop
  tgt_range.Begin();
  ref_range.Begin();
  while (tgt_range.NextInRange(&tgtRes)) {
    ref_range.NextInRange(&refRes);
    // Check if either the residue num or the reference residue num out of range.
    if ( tgtRes < 1 || tgtRes > NumResidues) {
      mprintf("    Warning: Rmsd: perres: Specified residue # %i is out of range.\n",tgtRes);
      continue;
    }
    if ( refRes < 1 || refRes > NumResidues ) {
      mprintf("    Warning: Rmsd: perres: Specified reference residue # %i is out of range.\n",
              refRes);
      continue;
    }
    N++;
    // Create dataset for res - if already present this returns NULL
    DataSet *prDataSet = PerResRMSD->AddMultiN(DataSet::DOUBLE, "", 
                                               currentParm->ResidueName(tgtRes-1),
                                               tgtRes);
    if (prDataSet != NULL) DFL->Add(perresout, prDataSet);

    // Setup mask strings. Note that masks are based off user residue nums
    sprintf(tgtArg,":%i%s",tgtRes,perresmask);
    tgtResMask[N].SetMaskString(tgtArg);
    sprintf(refArg,":%i%s",refRes,perresmask);
    refResMask[N].SetMaskString(refArg);
    //mprintf("DEBUG: RMSD: PerRes: Mask %s RefMask %s\n",tgtArg,refArg);

    // Setup the reference mask
    if (RefParm->SetupIntegerMask(refResMask[N])) {
      mprintf("      perres: Could not setup reference mask for residue %i\n",refRes);
      continue;
    }
    if (refResMask[N].None()) {
      mprintf("      perres: No atoms selected for reference residue %i\n",refRes);
      continue;
    }

    // Setup the target mask
    if (currentParm->SetupIntegerMask(tgtResMask[N])) {
      mprintf("      perres: Could not setup target mask for residue %i\n",tgtRes);
      continue;
    }
    if (tgtResMask[N].None()) {
      mprintf("      perres: No atoms selected for target residue %i\n",tgtRes);
      continue;
    }

    // Check that # atoms in target and reference masks match
    if (tgtResMask[N].Nselected() != refResMask[N].Nselected()) {
      mprintf("      perres: Res %i: # atoms in Tgt [%i] != # atoms in Ref [%i]\n",
              tgtRes,tgtResMask[N].Nselected(),refResMask[N].Nselected());
      continue;
    }

    // Indicate that these masks were properly set up
    resIsActive[N]=true;
  }   

  // Check pointer to the output file
  if (perresout!=NULL) {
    if (DFL->GetDataFile(perresout)==NULL) {
      mprinterr("Error: RMSD: Perres output file could not be set up.\n");
      return 1;
    }
  }

  // Allocate memory for residue frame and residue reference frame. The size 
  // of each Frame is initially allocated to the maximum number of atoms.
  // Although initial masses are wrong this is ok since the number of atoms 
  // and masses will change when residue RMSD is actually being calcd.
  if (ResRefFrame!=NULL) delete ResRefFrame;
  ResRefFrame = new Frame( RefParm->FindResidueMaxNatom(), RefParm->Mass() );
  //ResRefFrame->Info("ResRefFrame");
  if (ResFrame!=NULL) delete ResFrame;
  ResFrame = new Frame( currentParm->FindResidueMaxNatom(), currentParm->Mass() );
  //ResFrame->Info("ResFrame");

  return 0;
}

// Rmsd::setup()
/** Called every time the trajectory changes. Set up FrameMask for the new 
  * parmtop and allocate space for selected atoms from the Frame.
  */
int Rmsd::setup() {

  if ( currentParm->SetupIntegerMask( FrameMask ) ) return 1;
  if ( FrameMask.None() ) {
    mprintf("    Error: Rmsd::setup: No atoms in mask.\n");
    return 1;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  SelectedFrame.SetupFrameFromMask(FrameMask, currentParm->Mass());

  // Reference setup
  if (RefSetup( currentParm )) return 1;
  
  // Check that num atoms in frame mask from this parm match ref parm mask
  if ( RefNselected() != FrameMask.Nselected() ) {
    mprintf( "    Warning: Number of atoms in RMS mask (%i) does not equal number of\n",
            FrameMask.Nselected());
    mprintf( "             atoms in reference mask (%i).\n",RefNselected());
    return 1;
  }

  // Per residue rmsd setup
  if (perres) { 
    if (this->perResSetup(GetRefParm())) return 1;
  }

  if (!isSeparate_)
    mprintf("\t%i atoms selected.\n",FrameMask.Nselected());

  return 0;
}

// Rmsd::action()
/** Called every time a frame is read in. Calc RMSD. If not first and not
  * RefTraj SetRefStructure has already been called. When fitting, 
  * SetRefStructure pre-centers the reference coordinates at the origin
  * and puts the translation from origin to reference in Trans[3-5]. 
  */
int Rmsd::action() {
  double R, U[9];

  // Perform any needed reference actions
  RefAction(currentFrame, Trans+3);

  // Set selected frame atoms. Masses have already been set.
  SelectedFrame.SetCoordinates(*currentFrame, FrameMask);

  // DEBUG
/*  mprintf("  DEBUG: RMSD: First atom coord in SelectedFrame is : "); 
  SelectedFrame.printAtomCoord(0);
  mprintf("  DEBUG: RMSD: First atom coord in SelectedRef is : ");
  SelectedRef.printAtomCoord(0);
*/

  if (nofit) {
    R = SelectedFrame.RMSD(&SelectedRef_, useMass_);
  } else {
    R = SelectedFrame.RMSD_CenteredRef(SelectedRef_, U, Trans, useMass_);
    currentFrame->Trans_Rot_Trans(Trans,U);
  }

  rmsd->Add(frameNum, &R);

  // ---=== Per Residue RMSD ===---
  // Set reference and selected frame for each residue using the previously
  // set-up masks in refResMask and tgtResMask. Use SetFrameFromMask instead
  // of SetFrameCoordsFromMask since each residue can be a different size.
  if (perres) {
    for (int N=0; N < NumResidues; N++) {
      if (!resIsActive[N]) {
        //mprintf("DEBUG:           [%4i] Not Active.\n",N);
        continue;
      }
      ResRefFrame->SetFrame(RefFrame_, refResMask[N]);
      ResFrame->SetFrame(*currentFrame, tgtResMask[N]);
      if (perrescenter) {
        ResFrame->ShiftToGeometricCenter( );
        ResRefFrame->ShiftToGeometricCenter( );
      }
      R = ResFrame->RMSD(ResRefFrame,useMass_);
      //mprintf("DEBUG:           [%4i] Res [%s] nofit RMSD to [%s] = %lf\n",N,
      //        tgtResMask[N]->MaskString(),refResMask[N]->MaskString(),R);
      // NOTE: Should check for error on AddData?
      PerResRMSD->AddData(frameNum, &R, N);
    }
  }

  return 0;
}

// Rmsd::print()
/** For per-residue RMSD only. Sync the per-residue RMSD data set since
  * it is not part of the master DataSetList in Cpptraj. Setup output
  * file options. Calculate averages if requested.
  */
void Rmsd::print() {
  DataFile *outFile;

  if (!perres || PerResRMSD==NULL) return;
  // Sync dataset list here since it is not part of master dataset list
  PerResRMSD->Sync();
  // Per-residue output
  outFile = DFL->GetDataFile(perresout);
  if (outFile!=NULL) {
    // Set output file to be inverted if requested
    if (perresinvert) 
      outFile->ProcessArgs("invert");
    mprintf("    RMSD: Per-residue: Writing data for %i residues to %s\n",
            PerResRMSD->Size(), outFile->Filename());
  }

  // Average
  if (perresavg==NULL) return;
  int Nperres = PerResRMSD->Size();
  // Use the per residue rmsd dataset list to add one more for averaging
  DataSet *PerResAvg = PerResRMSD->Add(DataSet::DOUBLE, (char*)"AvgRMSD", "AvgRMSD");
  // another for stdev
  DataSet *PerResStdev = PerResRMSD->Add(DataSet::DOUBLE, (char*)"Stdev", "Stdev");
  // Add the average and stdev datasets to the master datafile list
  outFile = DFL->Add(perresavg, PerResAvg);
  outFile = DFL->Add(perresavg, PerResStdev);
  outFile->ProcessArgs("xlabel Residue");
  // For each residue, get the average rmsd
  double stdev = 0;
  double avg = 0;
  for (int pridx = 0; pridx < Nperres; pridx++) {
    DataSet *tempDS = PerResRMSD->GetDataSetN(pridx);
    avg = tempDS->Avg(&stdev);
    int dsidx = tempDS->Idx() - 1; // When set up actual resnum is used - change?
    PerResAvg->Add(dsidx, &avg);
    PerResStdev->Add(dsidx,&stdev);
  }
}
 

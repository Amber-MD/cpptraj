// RMSD
#include <cstdio> // for sprintf
#include "Action_Rmsd.h"
#include "CpptrajStdio.h"
#include "DataSet_double.h" // for SeparateInit

// CONSTRUCTOR
Rmsd::Rmsd() {
  rmsd=NULL;
  PerResRMSD=NULL;
  NumResidues=0;
  perresout=NULL;
  first=false;
  nofit=false;
  useMass=false;
  perres=false;
  RefTraj=NULL;
  RefParm=NULL;
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
  // If RMS to reference traj, close the traj. 
  if (RefTraj!=NULL) {
    RefTraj->EndTraj();
    delete RefTraj;
  }
  if (ResFrame!=NULL) delete ResFrame;
  if (ResRefFrame!=NULL) delete ResRefFrame;
  if (PerResRMSD!=NULL) delete PerResRMSD;
  // If separate, clean up the dataset
  if (isSeparate) delete rmsd;
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
    
// Rmsd::SetRefMask()
/** Setup reference mask based on maskRef. Requires RefParm to be set. Should 
  * only be called once.
  * If reference, this is called from init. If first, this is called from setup.
  */
int Rmsd::SetRefMask() {
  if (RefParm->SetupIntegerMask( RefMask, activeReference )) return 1;
  if (RefMask.None()) {
    mprintf("    Error: Rmsd::SetRefMask: No atoms in reference mask.\n");
    return 1;
  }
  // Check if reference parm has masses
  if (useMass && RefParm->mass==NULL) {
    mprintf("    Warning: usemass: Ref Parmtop %s does not contain mass info.\n",
            RefParm->parmName);
    mprintf("             Geometric center will be used instead.\n");
    useMass=false;
  }
  // Allocate frame for selected reference atoms
  SelectedRef.SetupFrameFromMask(&RefMask, RefParm->mass);
  //mprintf("DEBUG: RefMask has %i atoms\n",RefMask.Nselected);
  return 0;
}

// Rmsd::SetRefStructure()
/** Set coordinates of SelectedRef according to RefMask. If performing 
  * fitting, pre-translate the reference to the origin and save the
  * translation from origin to original reference location in
  * Trans[3-5]. When first, this is called from action once. When
  * RefTraj, this is called from action each time a reference structure
  * is read from RefTraj. Otherwise called from init once.
  */
void Rmsd::SetRefStructure() {
  SelectedRef.SetFrameCoordsFromMask(RefFrame.X, &RefMask);
  if (!nofit) 
    SelectedRef.CenterReference(Trans+3, useMass);
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
  char *referenceName, *mask0, *maskRef, *reftraj;
  char *rmsdFile;
  int refindex, referenceKeyword;
  Frame *TempFrame = NULL;

  // Check for keywords
  referenceKeyword=actionArgs.hasKey("reference"); // For compatibility with ptraj
  referenceName=actionArgs.getKeyString("ref",NULL);
  refindex=actionArgs.getKeyInt("refindex",-1);
  reftraj = actionArgs.getKeyString("reftraj",NULL);
  if (reftraj!=NULL) {
    RefParm = PFL->GetParm(actionArgs);
    if (RefParm==NULL) {
      mprinterr("Error: Rmsd: Could not get parm for reftraj %s.\n",reftraj);
      return 1;
    }
  }
  nofit = actionArgs.hasKey("nofit");
  first = actionArgs.hasKey("first");
  useMass = actionArgs.hasKey("mass");
  rmsdFile = actionArgs.getKeyString("out",NULL);
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
  mask0 = actionArgs.getNextMask();
  FrameMask.SetMaskString(mask0);
  // Get RMS mask string for reference
  maskRef = actionArgs.getNextMask();
  // If no reference mask specified, make same as RMS mask
  if (maskRef==NULL) maskRef=mask0; 
  RefMask.SetMaskString(maskRef);

  // Set up the RMSD data set. 
  rmsd = DSL->Add(DOUBLE, actionArgs.getNextString(),"RMSD");
  if (rmsd==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(rmsdFile,rmsd);

  // Check reference structure
  if (!first && referenceName==NULL && refindex==-1 && referenceKeyword==0 && reftraj==NULL) {
    mprintf("    Warning: Rmsd::init: No reference structure given. Defaulting to first.\n");
    first=true;
  }

  // If not using first frame, set up reference now.
  if (!first) {
    // Check if reference will be a series of frames from a trajectory
    if (reftraj!=NULL) {
      // Attempt to set up reference trajectory
      RefTraj = new TrajectoryFile();
      if (RefTraj->SetupRead(reftraj, NULL, RefParm)) {
        mprinterr("Error: Rmsd: Could not set up reftraj %s.\n",reftraj);
        delete RefTraj;
        RefTraj=NULL;
        return 1;
      } 
      RefFrame.SetupFrameV(RefParm->natom, RefParm->mass, RefTraj->HasVelocity());
      // Set up reference mask. Ref structures will be read in during action
      if ( SetRefMask() ) return 1; 
    } else {
      // Attempt to get reference index by name/tag
      if (referenceName!=NULL)
        refindex=FL->GetFrameIndex(referenceName);

      // For compatibility with ptraj, if 'reference' specified use first 
      // specified reference.
      if (referenceKeyword) refindex=0;

      // Get reference frame by index
      TempFrame=FL->GetFrame(refindex);
      if (TempFrame==NULL) {
        mprinterr("    Error: Rmsd::init: Could not get reference index %i\n",refindex);
        return 1;
      }
      RefFrame = *TempFrame;
      // Set reference parm
      RefParm=FL->GetFrameParm(refindex);
      // Set up reference mask and structure
      if ( SetRefMask() ) return 1;
      SetRefStructure();
    }
    //RefFrame.printAtomCoord(0);
    //fprintf(stderr,"  NATOMS IN REF IS %i\n",RefFrame.natom); // DEBUG
  }

  //rmsd->Info();
  mprintf("    RMSD: (%s), reference is ",FrameMask.MaskString());
  if (reftraj!=NULL) {
    // Set up reference trajectory and open
    mprintf("trajectory %s with %i frames",RefTraj->TrajName(),RefTraj->Total_Read_Frames());
    if (RefTraj->BeginTraj(false)) {
      mprinterr("Error: Rmsd: Could not open reference trajectory.\n");
      return 1;
    }
  } else if (first)
    mprintf("first frame");
  else if (referenceName!=NULL)
    mprintf("%s",referenceName);
  else
    mprintf("reference index %i",refindex);
  mprintf(" (%s)",RefMask.MaskString());
  if (nofit)
    mprintf(", no fitting");
  else
    mprintf(", with fitting");
  if (useMass) 
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
  isSeparate = true;
  debug = debugIn;
  useMass = massIn;
  // Also set useMassOriginalValue since this is NOT called from 
  // Init.
  useMassOriginalValue = useMass;
  // Only first for reference for now
  first = true;
  RefParm = NULL;

  // Set the RMS mask string for target and reference
  FrameMask.SetMaskString(mask0);
  RefMask.SetMaskString(mask0);
  //mprintf("DEBUG: Mask [%s] RefMask [%s]\n",FrameMask.MaskString(),RefMask.MaskString());

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
int Rmsd::perResSetup() {
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
    // Setup dataset name for this residue
    currentParm->ResName(tgtArg,tgtRes-1);
    // Create dataset for res - if already present this returns NULL
    DataSet *prDataSet = PerResRMSD->AddMultiN(DOUBLE, "", currentParm->ResidueName(tgtRes-1),
                                               tgtRes);
    if (prDataSet != NULL) DFL->Add(perresout, prDataSet);

    // Setup mask strings. Note that masks are based off user residue nums
    sprintf(tgtArg,":%i%s",tgtRes,perresmask);
    tgtResMask[N].SetMaskString(tgtArg);
    sprintf(refArg,":%i%s",refRes,perresmask);
    refResMask[N].SetMaskString(refArg);
    //mprintf("DEBUG: RMSD: PerRes: Mask %s RefMask %s\n",tgtArg,refArg);

    // Setup the reference mask
    if (RefParm->SetupIntegerMask(refResMask[N], activeReference )) {
      mprintf("      perres: Could not setup reference mask for residue %i\n",refRes);
      continue;
    }
    if (refResMask[N].None()) {
      mprintf("      perres: No atoms selected for reference residue %i\n",refRes);
      continue;
    }

    // Setup the target mask
    if (currentParm->SetupIntegerMask(tgtResMask[N], activeReference)) {
      mprintf("      perres: Could not setup target mask for residue %i\n",tgtRes);
      continue;
    }
    if (tgtResMask[N].None()) {
      mprintf("      perres: No atoms selected for target residue %i\n",tgtRes);
      continue;
    }

    // Check that # atoms in target and reference masks match
    if (tgtResMask[N].Nselected != refResMask[N].Nselected) {
      mprintf("      perres: Res %i: # atoms in Tgt [%i] != # atoms in Ref [%i]\n",
              tgtRes,tgtResMask[N].Nselected,refResMask[N].Nselected);
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
  // The number of atoms and masses will change based on which residue is 
  // currently being calcd.
  if (ResRefFrame!=NULL) delete ResRefFrame;
  ResRefFrame = new Frame();
  ResRefFrame->SetupFrame(RefParm->natom, RefParm->mass);
  if (ResFrame!=NULL) delete ResFrame;
  ResFrame = new Frame();
  ResFrame->SetupFrame(currentParm->natom, currentParm->mass);

  return 0;
}

// Rmsd::setup()
/** Called every time the trajectory changes. Set up FrameMask for the new 
  * parmtop and allocate space for selected atoms from the Frame.
  */
int Rmsd::setup() {

  if ( currentParm->SetupIntegerMask( FrameMask, activeReference ) ) return 1;
  if ( FrameMask.None() ) {
    mprintf("    Error: Rmsd::setup: No atoms in mask.\n");
    return 1;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  SelectedFrame.SetupFrameFromMask(&FrameMask, currentParm->mass);

  // first: If RefParm not set, set it here and setup the reference mask
  if (first && RefParm==NULL) {
    RefParm = currentParm;
    if ( SetRefMask( ) ) return 1;
  }
  
  // Check that num atoms in frame mask from this parm match ref parm mask
  if ( RefMask.Nselected != FrameMask.Nselected ) {
    mprintf( "    Warning: Number of atoms in RMS mask (%i) does not equal number of\n",
            FrameMask.Nselected);
    mprintf( "             atoms in reference mask (%i).\n",RefMask.Nselected);
    return 1;
  }

  // Per residue rmsd setup
  if (perres) { 
    if (this->perResSetup()) return 1;
  }

  if (!isSeparate)
    mprintf("\t%i atoms selected.\n",FrameMask.Nselected);

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

  if (first) {
    RefFrame = *currentFrame;
    first = false;
    SetRefStructure();
  }

  // reftraj: Get the next frame from the reference trajectory
  //          If no more frames are left, the last frame will be used. This
  //          could eventually be changed so that the trajectory loops.
  if (RefTraj!=NULL) {
    //mprintf("DBG: RMSD reftraj: Getting ref traj frame %i\n",RefTraj->front()->CurrentFrame());
    // NOTE: If there are no more frames in the trajectory the frame should
    //       remain on the last read frame. Close and reopen? Change ref?
    RefTraj->GetNextFrame(RefFrame);
    SetRefStructure(); 
  }

  // Set selected frame atoms. Masses have already been set.
  SelectedFrame.SetFrameCoordsFromMask(currentFrame->X, &FrameMask);

  // DEBUG
/*  mprintf("  DEBUG: RMSD: First atom coord in SelectedFrame is : "); 
  SelectedFrame.printAtomCoord(0);
  mprintf("  DEBUG: RMSD: First atom coord in SelectedRef is : ");
  SelectedRef.printAtomCoord(0);
*/

  if (nofit) {
    R = SelectedFrame.RMSD(&SelectedRef, useMass);
  } else {
    //R = SelectedFrame.RMSD(&SelectedRef, U, Trans, useMass);
    R = SelectedFrame.RMSD_CenteredRef(SelectedRef, U, Trans, useMass);
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
      ResRefFrame->SetFrameFromMask(&RefFrame, &(refResMask[N]));
      ResFrame->SetFrameFromMask(currentFrame, &(tgtResMask[N]));
      if (perrescenter) {
        ResFrame->ShiftToGeometricCenter( );
        ResRefFrame->ShiftToGeometricCenter( );
      }
      R = ResFrame->RMSD(ResRefFrame,useMass);
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
  * it is not part of the master DataSetList in CpptrajState. Setup output
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
      outFile->SetInverted();
    mprintf("    RMSD: Per-residue: Writing data for %i residues to %s\n",
            PerResRMSD->Size(), outFile->Filename());
  }

  // Average
  if (perresavg==NULL) return;
  int Nperres = PerResRMSD->Size();
  // Use the per residue rmsd dataset list to add one more for averaging
  DataSet *PerResAvg = PerResRMSD->Add(DOUBLE, (char*)"AvgRMSD", "AvgRMSD");
  // another for stdev
  DataSet *PerResStdev = PerResRMSD->Add(DOUBLE, (char*)"Stdev", "Stdev");
  // Add the average and stdev datasets to the master datafile list
  outFile = DFL->Add(perresavg, PerResAvg);
  outFile = DFL->Add(perresavg, PerResStdev);
  outFile->SetXlabel((char*)"Residue");
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
 

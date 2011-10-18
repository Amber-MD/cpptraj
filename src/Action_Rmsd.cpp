// RMSD
#include <cstdio> // for sprintf
#include "Action_Rmsd.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Rmsd::Rmsd() {
  rmsd=NULL;
  PerResRMSD=NULL;
  nres=0;
  perresout=NULL;
  nofit=false;
  first=false;
  useMass=false;
  perres=false;
  RefTraj=NULL;
  RefFrame=NULL;
  RefParm=NULL;
  SelectedRef=NULL;
  SelectedFrame=NULL;
  ResFrame=NULL;
  ResRefFrame=NULL;
  perresmask=NULL;
  perrescenter=false;
  perresinvert=false;
}

// DESTRUCTOR
Rmsd::~Rmsd() {
  //mprinterr("RMSD DESTRUCTOR\n");
  // If first, ref Frame was allocd (not assigned from reference Frame List)
  if (first && RefFrame!=NULL)
    delete RefFrame;
  if (RefTraj!=NULL) {
    RefTraj->EndTraj();
    delete RefTraj;
    if (RefFrame!=NULL) delete RefFrame;
  }
  if (SelectedRef!=NULL) delete SelectedRef;
  if (SelectedFrame!=NULL) delete SelectedFrame;
  if (ResFrame!=NULL) delete ResFrame;
  if (ResRefFrame!=NULL) delete ResRefFrame;
  if (PerResRMSD!=NULL) delete PerResRMSD;
  // Free up perres masks
  std::vector<AtomMask*>::iterator mask;
  for (mask = tgtResMask.begin(); mask != tgtResMask.end(); mask++) 
    delete (*mask);
  for (mask = refResMask.begin(); mask != refResMask.end(); mask++) 
    delete (*mask);
}

/* Rmsd::resizeResMasks()
 * For perres rmsd. If the current number of residues is greater than
 * the size of the residue mask lists, allocate as many extra masks
 * as needed. 
 */
void Rmsd::resizeResMasks() {
  int currentSize = (int)tgtResMask.size();
  if (nres > currentSize) {
    tgtResMask.resize(nres, (AtomMask*)NULL);
    refResMask.resize(nres, (AtomMask*)NULL);
    for (int res=currentSize; res < nres; res++) {
      tgtResMask[res]=new AtomMask();
      refResMask[res]=new AtomMask();
    }
  }
} 
    
/* Rmsd::SetRefMask()
 * Setup reference mask based on maskRef. Requires RefParm to be set. Should 
 * only be called once.
 * If reference, this is called from init. If first, this is called from setup.
 */
int Rmsd::SetRefMask() {
  if ( RefMask.SetupMask(RefParm,debug) ) return 1;
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
  SelectedRef = new Frame(&RefMask, RefParm->mass);

  return 0;
}

/* Rmsd::init()
 * Called once before traj processing. Set up reference info.
 * Expected call: 
 * rmsd <name> <mask> [<refmask>] [out filename] [nofit] [mass]
 *      [ first | ref <filename> | refindex <#> | 
          reftraj <filename> [parm <parmname> | parmindex <#>] ] 
 *      [ perres perresout <filename> [range <res range>] [refrange <ref res range>] 
 *        [perresmask <addtl mask>] [perresinvert] [perrescenter] ]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Rmsd::init( ) {
  char *referenceName, *mask0, *maskRef, *reftraj;
  char *rmsdFile;
  int refindex, referenceKeyword;

  // Check for keywords
  referenceKeyword=A->hasKey("reference"); // For compatibility with ptraj
  referenceName=A->getKeyString("ref",NULL);
  refindex=A->getKeyInt("refindex",-1);
  reftraj = A->getKeyString("reftraj",NULL);
  if (reftraj!=NULL) {
    RefParm = PFL->GetParm(A);
    if (RefParm==NULL) {
      mprinterr("Error: Rmsd: Could not get parm for reftraj %s.\n",reftraj);
      return 1;
    }
  }
  nofit = A->hasKey("nofit");
  first = A->hasKey("first");
  useMass = A->hasKey("mass");
  perres = A->hasKey("perres");
  rmsdFile = A->getKeyString("out",NULL);
  // Per-res keywords
  perresout = A->getKeyString("perresout",NULL);
  perresinvert = A->hasKey("perresinvert");
  ResRange.SetRange( A->getKeyString("range",NULL) );
  RefRange.SetRange( A->getKeyString("refrange",NULL) );
  perresmask = A->getKeyString("perresmask",(char*)"");
  perrescenter = A->hasKey("perrescenter");

  // Get the RMS mask string for frames
  mask0 = A->getNextMask();
  FrameMask.SetMaskString(mask0);
  // Get RMS mask string for reference
  maskRef = A->getNextMask();
  // If no reference mask specified, make same as RMS mask
  if (maskRef==NULL) maskRef=mask0; 
  RefMask.SetMaskString(maskRef);

  // Set up the RMSD data set
  rmsd = DSL->Add(DOUBLE, A->getNextString(),"RMSD");
  if (rmsd==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(rmsdFile,rmsd);

  if (!first && referenceName==NULL && refindex==-1 && referenceKeyword==0 && reftraj==NULL) {
    mprintf("    Warning: Rmsd::init: No reference structure given. Defaulting to first.\n");
    first=true;
  }

  if (!first) {
    // Check if reference will be a series of frames from a trajectory
    if (reftraj!=NULL) {
      if ( SetRefMask() ) return 1;
      // Attempt to set up reference trajectory
      RefTraj = new TrajectoryFile();
      if (RefTraj->SetupRead(reftraj, NULL, RefParm)) {
        mprinterr("Error: Rmsd: Could not set up reftraj %s.\n",reftraj);
        delete RefTraj;
        RefTraj=NULL;
        return 1;
      } 
      RefFrame = new Frame(RefParm->natom, RefParm->mass, RefTraj->HasVelocity());
    } else {
      // Attempt to get reference index by name
      if (referenceName!=NULL)
        refindex=FL->GetFrameIndex(referenceName);

      // For compatibility with ptraj, if 'reference' specified use first reference structure
      if (referenceKeyword) refindex=0;

      // Get reference frame by index
      RefFrame=FL->GetFrame(refindex);
      if (RefFrame==NULL) {
        mprinterr("    Error: Rmsd::init: Could not get reference index %i\n",refindex);
        return 1;
      }
      // Set reference parm
      RefParm=FL->GetFrameParm(refindex);
      // Setup reference mask here since reference frame/parm are allocated
      if ( SetRefMask() ) return 1;
      //RefFrame->printAtomCoord(0);
      //fprintf(stderr,"  NATOMS IN REF IS %i\n",RefFrame->P->natom); // DEBUG
    }
  }

  //rmsd->Info();
  mprintf("    RMSD: (%s), reference is ",FrameMask.maskString);
  if (reftraj!=NULL) {
    // Set up reference trajectory and open
    mprintf("trajectory %s with %i frames",RefTraj->TrajName(),RefTraj->Total_Read_Frames());
    if (RefTraj->BeginTraj(false)) {
      mprinterr("Error: Rmsd: Could not open reference trajectory.\n");
      return 1;
    }
  } else if (RefFrame==NULL)
    mprintf("first frame");
  else if (referenceName!=NULL)
    mprintf("%s",referenceName);
  else
    mprintf("reference index %i",refindex);
  mprintf(" (%s)",RefMask.maskString);
  if (nofit)
    mprintf(", no fitting");
  else
    mprintf(", with fitting");
  mprintf(".\n");
  if (useMass) 
    mprintf("          Center of mass will be used.\n");
  else
    mprintf("          Geometric center will be used.\n");

  // Per-residue RMSD
  if (perres) {
    mprintf("          No-fit RMSD will also be calculated for ");
    if (ResRange.Empty()) 
      mprintf("each solute residue");
    else
      mprintf("residues %s",ResRange.RangeArg());
    if (!RefRange.Empty())
      mprintf(" (reference residues %s)",RefRange.RangeArg());
    mprintf(" using mask [:X%s].\n",perresmask);
    if (perresout==NULL) {
      mprintf("Error: perres specified but no output filename given (perresout).\n");
      perres=false;
      return 1;
    }
    mprintf("          Per-residue output file is %s\n",perresout);
    if (perrescenter)
      mprintf("          perrescenter: Each residue will be centered prior to RMS calc.\n");
    if (perresinvert)
      mprintf("          perresinvert: Frames will be written in rows instead of columns.\n");
  }

  return 0;
}

/* Rmsd::perResSetup()
 * Perform setup required for per residue rmsd calculation.
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

  // If no target range previously specified do all solute residues
  if (ResRange.Empty()) {
    if (P->finalSoluteRes>0)
      nres = P->finalSoluteRes;
    else
      nres = P->nres; 
    tgt_range.SetRange(1,nres+1);
  } else
    tgt_range.SetRange(&ResRange);

  // If the reference range is empty, set it to match the target range
  if (RefRange.Empty()) 
    ref_range.SetRange(&tgt_range);
  else
    ref_range.SetRange(&RefRange);

  // Check that the number of reference residues matches number of target residues
  nres = tgt_range.Size();
  if (nres != ref_range.Size()) {
    mprinterr("Error: RMSD: PerRes: Number of residues %i does not match\n",nres);
    mprinterr("       number of reference residues %i.\n",ref_range.Size());
    return 1;
  }

  // Setup a dataset, target mask, and reference mask, for each residue.
  // Since we will only calculate per res rmsd for residues that can be
  // successfully set up, keep track of that as well.
  //mprinterr("DEBUG: Setting up %i masks and data for %s\n",nres,P->parmName);
  resizeResMasks();
  if (PerResRMSD==NULL) PerResRMSD=new DataSetList();
  resIsActive.reserve(nres);
  resIsActive.assign(nres,false);
  N = -1; // Set to -1 since increment is at top of loop
  tgt_range.Begin();
  ref_range.Begin();
  while (tgt_range.NextInRange(&tgtRes)) {
    ref_range.NextInRange(&refRes);
    // Check if either the residue num or the reference residue num out of range.
    if ( tgtRes < 1 || tgtRes > P->nres) {
      mprintf("    Warning: Rmsd: perres: Specified residue # %i is out of range.\n",tgtRes);
      continue;
    }
    if ( refRes < 1 || refRes > P->nres ) {
      mprintf("    Warning: Rmsd: perres: Specified reference residue # %i is out of range.\n",
              refRes);
      continue;
    }
    N++;
    // Setup dataset name for this residue
    P->ResName(tgtArg,tgtRes-1);
    // Create dataset for res - if already present this returns NULL
    DataSet *prDataSet = PerResRMSD->AddIdx(DOUBLE, tgtArg, tgtRes);
    if (prDataSet != NULL) DFL->Add(perresout, prDataSet);

    // Setup mask strings. Note that masks are based off user residue nums
    sprintf(tgtArg,":%i%s",tgtRes,perresmask);
    tgtResMask[N]->SetMaskString(tgtArg);
    sprintf(refArg,":%i%s",refRes,perresmask);
    refResMask[N]->SetMaskString(refArg);
    //mprintf("DEBUG: RMSD: PerRes: Mask %s RefMask %s\n",tgtArg,refArg);

    // Setup the reference mask
    if (refResMask[N]->SetupMask(RefParm, debug)) {
      mprintf("      perres: Could not setup reference mask for residue %i\n",refRes);
      continue;
    }
    if (refResMask[N]->None()) {
      mprintf("      perres: No atoms selected for reference residue %i\n",refRes);
      continue;
    }

    // Setup the target mask
    if (tgtResMask[N]->SetupMask(P, debug)) {
      mprintf("      perres: Could not setup target mask for residue %i\n",tgtRes);
      continue;
    }
    if (tgtResMask[N]->None()) {
      mprintf("      perres: No atoms selected for target residue %i\n",tgtRes);
      continue;
    }

    // Check that # atoms in target and reference masks match
    if (tgtResMask[N]->Nselected != refResMask[N]->Nselected) {
      mprintf("      perres: Res %i: # atoms in Tgt [%i] != # atoms in Ref [%i]\n",
              tgtRes,tgtResMask[N]->Nselected,refResMask[N]->Nselected);
      continue;
    }

    // Indicate that these masks were properly set up
    resIsActive[N]=true;
  }   

  // Check pointer to the output file
  if (DFL->GetDataFile(perresout)==NULL) {
    mprintf("Error: RMSD: Perres output file could not be set up.\n");
    return 1;
  }

  // Allocate memory for residue frame and residue reference frame. The size 
  // of each Frame is initially allocated to the maximum number of atoms.
  // The number of atoms and masses will change based on which residue is 
  // currently being calcd.
  if (ResRefFrame!=NULL) delete ResRefFrame;
  ResRefFrame = new Frame(RefParm->natom, RefParm->mass);
  if (ResFrame!=NULL) delete ResFrame;
  ResFrame = new Frame(P->natom, P->mass);

  return 0;
}

/* Rmsd::setup()
 * Called every time the trajectory changes. Set up FrameMask for the new 
 *parmtop and allocate space for selected atoms from the Frame.
 */
int Rmsd::setup() {

  if ( FrameMask.SetupMask(P,debug) ) return 1;
  if ( FrameMask.None() ) {
    mprintf("    Error: Rmsd::setup: No atoms in mask.\n");
    return 1;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  if (SelectedFrame!=NULL) delete SelectedFrame;
  SelectedFrame = new Frame(&FrameMask, P->mass);
  
  // first: If RefParm not set, set it here and set the reference mask.
  //        Should only occur once.
  if (first && RefParm==NULL) {
    RefParm=P;
    if ( SetRefMask() ) return 1;
  } 

  // Check that num atoms in frame mask from this parm match ref parm mask
  if ( RefMask.Nselected != FrameMask.Nselected ) {
    mprintf( "    Error: Number of atoms in RMS mask (%i) does not \n",FrameMask.Nselected);
    mprintf( "           equal number of atoms in Ref mask (%i).\n",RefMask.Nselected);
    return 1;
  }

  // Per residue rmsd setup
  if (perres) { 
    if (this->perResSetup()) return 1;
  }

  return 0;
}

/* Rmsd::action()
 * Called every time a frame is read in. Calc RMSD. If RefFrame is NULL
 * at this point we want to set the first frame read in as reference.
 */
int Rmsd::action() {
  double R, U[9], Trans[6];

  // first: If Ref is NULL, allocate this frame as reference
  //        Should only occur once.
  // NOTE: For MPI this will currently result in different references between threads.
  if (first && RefFrame==NULL) 
    RefFrame = F->Copy();

  // reftraj: Get the next frame from the reference trajectory
  //          If no more frames are left, the last frame will be used. This
  //          could eventually be changed so that the trajectory loops.
  if (RefTraj!=NULL) {
    //mprintf("DBG: RMSD reftraj: Getting ref traj frame %i\n",RefTraj->front()->CurrentFrame());
    // NOTE: If there are no more frames in the trajectory the frame should
    //       remain on the last read frame. Close and reopen? Change ref?
    RefTraj->GetNextFrame(RefFrame->X, RefFrame->V, RefFrame->box, (&RefFrame->T)); 
  }

  // Set selected reference atoms - always done since RMS fit modifies SelectedRef 
  SelectedRef->SetFrameCoordsFromMask(RefFrame->X, &RefMask);

  // Set selected frame atoms. Masses have already been set.
  SelectedFrame->SetFrameCoordsFromMask(F->X, &FrameMask);

  // DEBUG
/*  mprintf("  DEBUG: RMSD: First atom coord in SelectedFrame is : "); 
  SelectedFrame->printAtomCoord(0);
  mprintf("  DEBUG: RMSD: First atom coord in SelectedRef is : ");
  SelectedRef->printAtomCoord(0);
*/

  if (nofit) {
    R = SelectedFrame->RMSD(SelectedRef, useMass);
  } else {
    R = SelectedFrame->RMSD(SelectedRef, U, Trans, useMass);
    F->Translate(Trans);
    F->Rotate(U);
    F->Translate(Trans+3);
  }

  rmsd->Add(currentFrame, &R);

  // ---=== Per Residue RMSD ===---
  // Set reference and selected frame for each residue using the previously
  // set-up masks in refResMask and tgtResMask. Use SetFrameFromMask instead
  // of SetFrameCoordsFromMask since each residue can be a different size.
  if (perres) {
    for (int N=0; N < nres; N++) {
      if (!resIsActive[N]) {
        //mprintf("DEBUG:           [%4i] Not Active.\n",N);
        continue;
      }
      ResRefFrame->SetFrameFromMask(RefFrame, refResMask[N]);
      ResFrame->SetFrameFromMask(F, tgtResMask[N]);
      if (perrescenter)
        ResFrame->ShiftToCenter(ResRefFrame);
      R = ResFrame->RMSD(ResRefFrame,useMass);
      //mprintf("DEBUG:           [%4i] Res [%s] nofit RMSD to [%s] = %lf\n",N,
      //        tgtResMask[N]->maskString,refResMask[N]->maskString,R);
      // NOTE: Should check for error on AddData?
      PerResRMSD->AddData(currentFrame, &R, N);
    }
  }

  return 0;
}

/* Rmsd::print()
 * For per-residue RMSD only. Sync the per-residue RMSD data set since
 * it is not part of the master DataSetList in CpptrajState. Setup output
 * file options.
 */
void Rmsd::print() {
  DataFile *outFile;

  if (!perres) return;
  outFile = DFL->GetDataFile(perresout);
  if (outFile==NULL) {
    mprinterr("Error: RMSD: PerRes: Could not get perresout file %s.\n",perresout);
    return;
  }
  if (PerResRMSD==NULL) return;
  // Sync dataset list here since it is not part of master dataset list
  PerResRMSD->Sync();
  // Set output file to be inverted if requested
  if (perresinvert) 
    outFile->SetInverted();

  mprintf("    RMSD: Per-residue: Writing data for %i residues to %s\n",
          PerResRMSD->Size(), outFile->filename);
}
 

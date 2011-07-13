#include <cstdio> // for sprintf
#include "Action_Rmsd.h"
#include "CpptrajStdio.h"
// RMSD
using namespace std;

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
  outFile=NULL;
}

// DESTRUCTOR
Rmsd::~Rmsd() {
  // If first, ref Frame was allocd (not assigned from reference Frame List)
  if (first && RefFrame!=NULL)
    delete RefFrame;
  if (RefTraj!=NULL) {
    RefTraj->EndTraj();
    delete RefTraj;
  }
  if (SelectedRef!=NULL) delete SelectedRef;
  if (SelectedFrame!=NULL) delete SelectedFrame;
  if (ResFrame!=NULL) delete ResFrame;
  if (ResRefFrame!=NULL) delete ResRefFrame;
  if (PerResRMSD!=NULL) delete PerResRMSD;
  clearPerResMask();
}

/* Rmsd::clearPerResMask()
 * Free the memory used by the atom masks in PerResMask.
 */
void Rmsd::clearPerResMask() {
  while ( !PerResMask.empty() ) {
    AtomMask *ResMask = PerResMask.back();
    delete ResMask;
    PerResMask.pop_back();
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
 *      [perres] [perresout <filename>] [range <res range>] [refrange <ref res range>] 
 *        [perresmask <addtl mask>] [perresinvert] [perrescenter]
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
    } else {
      // Attempt to get reference index by name
      if (referenceName!=NULL)
        refindex=FL->GetFrameIndex(referenceName);

      // For compatibility with ptraj, if 'reference' specified use first reference structure
      if (referenceKeyword) refindex=0;

      // Get reference frame by index
      RefFrame=FL->GetFrame(refindex);
      if (RefFrame==NULL) {
        mprintf("    Error: Rmsd::init: Could not get reference index %i\n",refindex);
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
    mprintf("trajectory %s with %i frames.\n",RefTraj->TrajName(),RefTraj->Total_Read_Frames());
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
    if (ResRange.empty()) 
      mprintf("each solute residue");
    else
      mprintf("residues %s",ResRange.RangeArg());
    if (!RefRange.empty())
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

/* Rmsd::setup()
 * Called every time the trajectory changes. Set up FrameMask for the new parmtop
 * and allocate space for selected atoms from the Frame.
 */
int Rmsd::setup() {
  char resArg[1024];
  char refArg[1024];
  AtomMask *ResMask;
  int refRes;;

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

  // Check Mass
  if (useMass && P->mass==NULL) {
    mprintf("    Warning: Rmsd::setup: usemass: Parmtop %s does not contain mass info.\n",
            P->parmName);
    mprintf("             Geometric center will be used instead.\n");
    useMass=false;
  }

  // --------------------====  PER RESIDUE RMSD OPTION ====---------------------
  // If perres was specified, need a data set for each residue
  // NOTE THAT ALL RESIDUES FROM INPUT SHOULD BE SHIFTED BY -1
  if (perres) {
    if (PerResRMSD!=NULL) {
      // This is the second parm for which perres is being called for. 
      // Potentially problematic since there is no guarantee each residue
      // matches up in each parmtop. Just print a warning for now.
      mprintf("    Warning: RMSD: perres option in use for more than 1 prmtop.\n");
      mprintf("             Residue names are not guaranteed to match.\n");
    } else {
      PerResRMSD = new DataSetList();
    }

    // If no range previously specified do all solute residues.
    if (ResRange.empty()) {
      if (P->finalSoluteRes>0)
        nres = P->finalSoluteRes;
      else
        nres=P->nres;
      // Has to match user input where residue nums start from 1
      // Since ArgToRange returns up to and including, want 1-nres
      ResRange.SetRange(1,nres+1);
    } 
    mprintf("      RMSD: PerRes: Setting up for %i residues.\n",(int)ResRange.size());

    // If the reference range is empty, set it to match the residue range
    if (RefRange.empty()) RefRange.assign(ResRange.begin(), ResRange.end());

    // Check that the number of reference residues matches number of residues
    if (RefRange.size() != ResRange.size()) {
      mprinterr("Error: RMSD: PerRes: Number of residues %i does not match\n",(int)ResRange.size());
      mprinterr("       number of reference residues %i.\n",(int)RefRange.size());
      return 1;
    }

    // For each residue specified in the range, set up an atom mask for selected
    // and reference atoms, along with a data set.
    // PerResMask will hold both masks, the ref mask followed by selected mask.
    // NOTE: res = *it - 1; res is the internal resnum, *it the user resnum
    clearPerResMask();
    std::list<int>::iterator refit = RefRange.begin();
    for (std::list<int>::iterator it=ResRange.begin(); it!=ResRange.end(); it++) {
      // Get corresponding reference resnum
      refRes = *refit;
      refit++;

      // Setup mask strings - masks are based off user residue nums
      sprintf(resArg,":%i%s",*it,perresmask);
      sprintf(refArg,":%i%s",refRes,perresmask);
      //mprintf("DEBUG: RMSD: PerRes: Mask %s RefMask %s\n",resArg,refArg);

      // Set up reference mask for this residue.
      ResMask=new AtomMask();
      ResMask->SetMaskString(refArg);
      if ( ResMask->SetupMask(RefParm,0) ) {
        mprintf("Warning: RMSD: PerRes: Could not set up reference for residue %i\n",*it);
        delete ResMask;
        continue;
      }
      if (ResMask->None()) {
        mprintf("Warning: RMSD: PerRes: No atoms in reference for residue %i\n",*it);
        delete ResMask;
        continue;
      }
      //RefNselected = RefMask->Nselected;
      PerResMask.push_back(ResMask);

      // Set up mask for this residue.
      // If unable make sure reference mask is popped as well
      ResMask=new AtomMask();
      ResMask->SetMaskString(resArg);
      if ( ResMask->SetupMask(P,0) ) { // NOTE: Allow debug value in here?
        mprintf("Warning: RMSD: PerRes: Could not set up mask for residue %i\n",*it);
        delete ResMask;
        ResMask = PerResMask.back();
        delete ResMask;
        PerResMask.pop_back();
        continue;
      }
      if (ResMask->None()) {
        mprintf("Warning: RMSD: PerRes: No atoms in mask for residue %i\n",*it);
        delete ResMask;
        ResMask = PerResMask.back();
        delete ResMask;
        PerResMask.pop_back();
        continue;
      }

      // Check that number of atoms selected in parm is same as reference
      if ( (PerResMask.back())->Nselected != ResMask->Nselected) {
        mprintf("Warning: RMSD: PerRes: # atoms in mask for residue %i (%i) not equal\n",
                *it,ResMask->Nselected);
        mprintf("                       to # atoms in reference mask (%i).\n",
                (PerResMask.back())->Nselected);
        delete ResMask;
        ResMask = PerResMask.back();
        delete ResMask;
        PerResMask.pop_back();
        continue;
      }
      PerResMask.push_back(ResMask);

      // DEBUG
      //mprintf("PERRES_RMS: Mask %s RefMask %s\n",(PerResMask.back())->maskString,
      //        (PerResMask.back())->maskString);

      // Setup dataset name for this residue
      P->ResName(resArg,(*it)-1);
      // Create dataset for res - if already present this returns NULL
      DataSet *prDataSet = PerResRMSD->AddIdx(DOUBLE, resArg, *it);
      if (prDataSet != NULL) DFL->Add(perresout, prDataSet);
    } // END loop over residues in range

    // Set up pointer to the output file
    outFile = DFL->GetDataFile(perresout);
    if (outFile==NULL) {
      mprintf("Error: RMSD: Perres output file could not be set up.\n");
      return 1;
    }
    // Allocate memory for residue frame and residue reference frame. The size 
    // of each Frame is initially allocated to the maximum number of atoms.
    // The number of atoms and masses will change based on which residue is 
    // currently being calcd.
    if (ResRefFrame==NULL)
      ResRefFrame = new Frame(RefParm->natom, RefParm->mass);
    if (ResFrame==NULL)
      ResFrame = new Frame(P->natom, P->mass);
  }

  return 0;
}

/* Rmsd::action()
 * Called every time a frame is read in. Calc RMSD. If RefFrame is NULL
 * at this point we want to set the first frame read in as reference.
 */
int Rmsd::action() {
  double R, U[9], Trans[6];
  vector<AtomMask*>::iterator mask;

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

  // Per Residue RMSD - Set reference and selected frame for each mask in 
  // PerResMask. PerResMask contains the reference mask followed by the
  // selected mask.
  // Use SetFrameFromMask since each residue can be a different size
  if (perres) {
    PerResRMSD->Begin(); 
    for (mask = PerResMask.begin(); mask!=PerResMask.end(); mask++) {
      ResRefFrame->SetFrameFromMask(RefFrame, (*mask));
      mask++;
      ResFrame->SetFrameFromMask(F, (*mask));
      if (perrescenter) 
        ResFrame->ShiftToCenter(ResRefFrame); 
      R = ResFrame->RMSD(ResRefFrame,useMass);
      //mprintf("DEBUG:           Res %i nofit RMSD = %lf\n",res,R);
      // NOTE: Should check for error on AddData?
      PerResRMSD->AddData(currentFrame, &R);
    }
  }

  return 0;
}

/* Rmsd::print()
 * Write out per-residue RMSD
 */
void Rmsd::print() {
  if (!perres || outFile==NULL) return;
  if (PerResRMSD==NULL) return;
  // Sync dataset list here since it is not part of master dataset list
  PerResRMSD->Sync();
  // Set output file to be inverted if requested
  if (perresinvert) 
    outFile->SetInverted();

  mprintf("    RMSD: Per-residue: Writing data for %i residues to %s\n",
          (int)PerResMask.size() / 2, outFile->filename);
}
 

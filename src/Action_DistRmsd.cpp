// DISTRMSD
#include <cstdio> // for sprintf
#include "Action_DistRmsd.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DistRmsd::DistRmsd() {
  drmsd=NULL;
  RefTraj=NULL;
  RefParm=NULL;
  first = false;
}

// DESTRUCTOR
DistRmsd::~DistRmsd() {
  //mprinterr("RMSD DESTRUCTOR\n");
  if (RefTraj!=NULL) {
    RefTraj->EndTraj();
    delete RefTraj;
  }
}

// DistRmsd::SetRefMask()
/** Setup reference mask based on maskRef. Requires RefParm to be set. Should 
  * only be called once.
  * If reference, this is called from init. If first, this is called from setup.
  */
int DistRmsd::SetRefMask() {
  if ( RefParm->SetupIntegerMask(RefMask, activeReference) ) return 1;
  if (RefMask.None()) {
    mprintf("    Error: DistRmsd::SetRefMask: No atoms in reference mask.\n");
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

  return 0;
}

// DistRmsd::init()
/** Called once before traj processing. Set up reference info.
  * Expected call: 
  * drmsd <name> <mask> [<refmask>] [out filename] 
  *       [ first | ref <filename> | refindex <#> | 
  *         reftraj <filename> [parm <parmname> | parmindex <#>] ] 
  */
int DistRmsd::init( ) {
  char *referenceName, *mask0, *maskRef, *reftraj;
  char *rmsdFile;
  int refindex, referenceKeyword;
  Frame *TempFrame=NULL;

  // Check for keywords
  referenceKeyword=actionArgs.hasKey("reference"); // For compatibility with ptraj
  referenceName=actionArgs.getKeyString("ref",NULL);
  refindex=actionArgs.getKeyInt("refindex",-1);
  reftraj = actionArgs.getKeyString("reftraj",NULL);
  if (reftraj!=NULL) {
    RefParm = PFL->GetParm(actionArgs);
    if (RefParm==NULL) {
      mprinterr("Error: DistRmsd: Could not get parm for reftraj %s.\n",reftraj);
      return 1;
    }
  }
  first = actionArgs.hasKey("first");
  rmsdFile = actionArgs.getKeyString("out",NULL);

  // Get the RMS mask string for frames
  mask0 = actionArgs.getNextMask();
  TgtMask.SetMaskString(mask0);
  // Get RMS mask string for reference
  maskRef = actionArgs.getNextMask();
  // If no reference mask specified, make same as RMS mask
  if (maskRef==NULL) maskRef=mask0; 
  RefMask.SetMaskString(maskRef);

  // Set up the RMSD data set
  drmsd = DSL->Add(DOUBLE, actionArgs.getNextString(),"DRMSD");
  if (drmsd==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(rmsdFile,drmsd);

  // Set up reference structure
  if (!first && referenceName==NULL && refindex==-1 && referenceKeyword==0 && reftraj==NULL) {
    mprintf("    Warning: DistRmsd::init: No reference structure given. Defaulting to first.\n");
    first=true;
  }

  // If not using first frame, set up reference now. 
  if (!first) {
    // Check if reference will be a series of frames from a trajectory
    if (reftraj!=NULL) {
      // Attempt to set up reference trajectory
      RefTraj = new TrajectoryFile();
      if (RefTraj->SetupRead(reftraj, NULL, RefParm)) {
        mprinterr("Error: DistRmsd: Could not set up reftraj %s.\n",reftraj);
        delete RefTraj;
        RefTraj=NULL;
        return 1;
      } 
      RefFrame.SetupFrameV(RefParm->natom, RefParm->mass, RefTraj->HasVelocity());
    // Attempt to get reference index by name/tag
    } else {
      if (referenceName!=NULL)
        refindex=FL->GetFrameIndex(referenceName);

      // For compatibility with ptraj, if 'reference' specified use first 
      // specified reference.
      if (referenceKeyword) refindex=0;

      // Get reference frame by index
      TempFrame=FL->GetFrame(refindex);
      if (TempFrame==NULL) {
        mprintf("    Error: DistRmsd::init: Could not get reference index %i\n",refindex);
        return 1;
      }
      RefFrame = *TempFrame;
      // Set reference parm
      RefParm=FL->GetFrameParm(refindex);
    }
    // Setup reference mask 
    if ( SetRefMask( ) ) return 1;
    //RefFrame.printAtomCoord(0);
    //fprintf(stderr,"  NATOMS IN REF IS %i\n",RefFrame.natom); // DEBUG
  }

  mprintf("    DISTRMSD: (%s), reference is ",TgtMask.MaskString());
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
  mprintf(" (%s).\n",RefMask.MaskString());

  return 0;
}

// DistRmsd::setup()
/** Called every time the trajectory changes. Set up TgtMask for the new 
  * parmtop and allocate space for selected atoms from the Frame.
  */
int DistRmsd::setup() {

  if ( currentParm->SetupIntegerMask(TgtMask, activeReference) ) return 1;
  if ( TgtMask.None() ) {
    mprintf("    Error: DistRmsd::setup: No atoms in mask.\n");
    return 1;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  SelectedTgt.SetupFrameFromMask(&TgtMask, currentParm->mass);

  // first: If RefParm not set, set it here and setup the reference mask
  if (first && RefParm==NULL) {
    RefParm = currentParm;
    if ( SetRefMask( ) ) return 1;
  }
  
  // Check that num atoms in frame mask from this parm match ref parm mask
  if ( RefMask.Nselected != TgtMask.Nselected ) {
    mprintf( "    Error: Number of atoms in RMS mask (%i) does not \n",TgtMask.Nselected);
    mprintf( "           equal number of atoms in Ref mask (%i).\n",RefMask.Nselected);
    return 1;
  }

  return 0;
}

// DistRmsd::action()
/** Called every time a frame is read in. Calc distance RMSD.
  * If first is true, set the first frame read in as reference.
  */
int DistRmsd::action() {
  double DR;

  if (first) {
    RefFrame = *currentFrame;
    first = false;
  }

  // reftraj: Get the next frame from the reference trajectory
  //          If no more frames are left, the last frame will be used. This
  //          could eventually be changed so that the trajectory loops.
  if (RefTraj!=NULL) {
    //mprintf("DBG: RMSD reftraj: Getting ref traj frame %i\n",RefTraj->front()->CurrentFrame());
    // NOTE: If there are no more frames in the trajectory the frame should
    //       remain on the last read frame. Close and reopen? Change ref?
    RefTraj->GetNextFrame(RefFrame); 
  }

  // Set selected reference atoms - always done since RMS fit modifies SelectedRef 
  SelectedRef.SetFrameCoordsFromMask(RefFrame.X, &RefMask);

  // Set selected frame atoms. Masses have already been set.
  SelectedTgt.SetFrameCoordsFromMask(currentFrame->X, &TgtMask);

  // DEBUG
/*  mprintf("  DEBUG: RMSD: First atom coord in SelectedTgt is : "); 
  SelectedTgt->printAtomCoord(0);
  mprintf("  DEBUG: RMSD: First atom coord in SelectedRef is : ");
  SelectedRef->printAtomCoord(0);
*/

  DR = SelectedTgt.DISTRMSD( &SelectedRef );

  drmsd->Add(frameNum, &DR);

  return 0;
}


/* Action: Rms2d
 * Perform RMS calculation between each input frame and each other input 
 * frame, or each frame read in from a separate reference traj and each 
 * input frame. 
 * The actual calcuation is performed in the print function.
 */
#include "Action_Rms2d.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#include <cstdio> //sprintf

// CONSTRUCTOR
Rms2d::Rms2d() {
  //fprintf(stderr,"Rms2d Con\n");
  nofit=false;
  useMass=false;
  rmsdFile=NULL;
  RefTraj=NULL;
  RefParm=NULL;
} 

// DESTRUCTOR
Rms2d::~Rms2d() { 
  if (RefTraj!=NULL) {
    RefTraj->front()->End();
    delete RefTraj;
  }
}

/*
 * Rms2d::init()
 * Expected call: rms2d <mask> <refmask> rmsout <filename> [mass] [nofit] 
                  [reftraj <traj> [parm <parmname> | parmindex <#>]] 
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Rms2d::init() {
  char *mask0, *maskRef, *reftraj;
  int temp = 0;

  // Get keywords
  nofit = A->hasKey("nofit");
  useMass = A->hasKey("mass");
  rmsdFile = A->getKeyString("rmsout",NULL);
  reftraj = A->getKeyString("reftraj",NULL);
  if (reftraj!=NULL) {
    RefParm = PFL->GetParm(A);
    if (RefParm==NULL) {
      mprinterr("Error: Rms2d: Could not get parm for reftraj %s.\n",reftraj);
      return 1;
    }
  }
  // Require an output filename
  if (rmsdFile==NULL) {
    mprinterr("Error: Rms2d: No output filename specified; use 'rmsout' keyword.\n");
    return 1;
  }

  // Get the RMS mask string for frames
  mask0 = A->getNextMask();
  FrameMask.SetMaskString(mask0);
  // Get RMS mask string for reference
  maskRef = A->getNextMask();
  // If no reference mask specified, make same as RMS mask
  if (maskRef==NULL) maskRef=mask0;
  RefMask.SetMaskString(maskRef);

  // Check if reference will be a series of frames from a trajectory
  if (reftraj!=NULL) {
    // Attempt to set up reference trajectory
    RefTraj = new TrajinList();
    if (RefTraj->AddTrajin(reftraj, RefParm)) {
      mprinterr("Error: Rms2d: Could not set up reftraj %s.\n",reftraj);
      delete RefTraj;
      RefTraj=NULL;
      return 1;
    }
  }

  mprintf("    RMS2D: (%s) to (%s)",FrameMask.maskString,RefMask.maskString);
  if (reftraj!=NULL) {
    mprintf(" reference is trajectory:\n        ");
    RefTraj->SetupFrames();
    mprintf("          ");
    RefTraj->front()->Begin(&temp,0);
    mprintf("         ");
  }
  if (nofit)
    mprintf(" (no fitting)");
  if (useMass)
    mprintf(" (mass-weighted)");
  if (rmsdFile!=NULL) 
    mprintf(" output to %s",rmsdFile);
  mprintf("\n");

  return 0;
}

/*
 * Rms2d::setup()
 * Not important for Rms2d, initial pass is only for storing frames.
 */
int Rms2d::setup() {
  return 0;  
}

/*
 * Rms2d::action()
 * Store current frame as a reference frame.
 */
int Rms2d::action() {
  Frame *fCopy;

  fCopy = F->Copy();
  ReferenceFrames.Add(fCopy,P);
  
  return 0;
} 

/*
 * Rms2d::print()
 * Perform the rms calculation of each frame to each other frame.
 */
void Rms2d::print() {
  Frame *RefFrame;
  AmberParm *TgtParm;
  Frame *TgtFrame;
  Frame *SelectedRef=NULL;
  Frame *SelectedTgt=NULL;
  int lastrefpindex=-1;
  int lasttgtpindex=-1;
  double R, U[9], Trans[6];
  int current=0;
  int max=0;
  int totalref=0;
  int res=0; // Currently only used as dummy space for getnextframe
  DataSet *rmsdata;
  char setname[256];

  if (RefTraj==NULL) {
    totalref = ReferenceFrames.NumFrames();
    max = ReferenceFrames.NumFrames() * ReferenceFrames.NumFrames();
    mprintf("  RMS2D: Calculating RMSDs between each frame (%i total).\n  ",max);
  } else {
    totalref = RefTraj->front()->total_read_frames;
    max = totalref * ReferenceFrames.NumFrames();
    mprintf("  RMS2D: Calculating RMSDs between each input frame and each reference\n"); 
    mprintf("         trajectory %s frame (%i total).\n  ",
            RefTraj->front()->trajfilename, max);
  }

  // Set up progress Bar
  ProgressBar *progress = new ProgressBar(max);

  // LOOP OVER REFERENCE FRAMES
  for (int nref=0; nref < totalref; nref++) {
    progress->Update(current);
    if (RefTraj==NULL) 
      RefParm = ReferenceFrames.GetFrameParm( nref );
    // If the current ref parm not same as last ref parm, reset reference mask
    if (RefParm->pindex != lastrefpindex) {
      if ( RefMask.SetupMask(RefParm,debug) ) {
        mprinterr("Error: Rms2d: Could not set up reference mask for %s\n",RefParm->parmName);
        if (SelectedRef!=NULL) delete SelectedRef;
        if (SelectedTgt!=NULL) delete SelectedTgt;
        return;
      }
      if ( SelectedRef!=NULL ) delete SelectedRef;
      SelectedRef = new Frame(&RefMask, RefParm->mass);
      lastrefpindex = RefParm->pindex;
    }
    // Get the current reference frame
    if (RefTraj==NULL)
      RefFrame = ReferenceFrames.GetFrame( nref );
    else {
      if ( RefTraj->front()->NextFrame(&res) )
        RefFrame = RefTraj->front()->F;
    }
    // Set up dataset for this reference frame
    sprintf(setname,"Frame_%i",nref+1);
    rmsdata = RmsData.Add(DOUBLE, setname, "Rms2d");
    DFL->Add(rmsdFile,rmsdata);

    // LOOP OVER TARGET FRAMES
    for (int nframe=0; nframe < ReferenceFrames.NumFrames(); nframe++) {
      TgtParm = ReferenceFrames.GetFrameParm( nframe );
      // If the current frame parm not same as last frame parm, reset frame mask
      if (TgtParm->pindex != lasttgtpindex) {
        if ( FrameMask.SetupMask(TgtParm,debug) ) {
          mprinterr("Error: Rms2d: Could not set up target mask for %s\n",TgtParm->parmName);
          if (SelectedRef!=NULL) delete SelectedRef;
          if (SelectedTgt!=NULL) delete SelectedTgt;
          return;
        }
        // Check that num atoms in mask are the same
        if (FrameMask.Nselected != RefMask.Nselected) {
          mprinterr("Error: Rms2d: Num atoms in target mask (%i) != num atoms in ref mask (%i)\n",
                    FrameMask.Nselected, RefMask.Nselected);
          if (SelectedRef!=NULL) delete SelectedRef;
          if (SelectedTgt!=NULL) delete SelectedTgt;
          return;
        }
        if ( SelectedTgt!=NULL ) delete SelectedTgt;
        SelectedTgt = new Frame(&FrameMask, TgtParm->mass);
        lasttgtpindex = TgtParm->pindex;
      }
      // Get the current target frame
      TgtFrame = ReferenceFrames.GetFrame( nframe );

      // Set selected reference atoms - always done since RMS fit modifies SelectedRef
      SelectedRef->SetFrameCoordsFromMask(RefFrame->X, &RefMask);
      // Set selected target atoms
      SelectedTgt->SetFrameCoordsFromMask(TgtFrame->X, &FrameMask);

      // Perform RMS calculation
      if (nofit) {
        R = SelectedTgt->RMSD(SelectedRef, useMass);
      } else {
        R = SelectedTgt->RMSD(SelectedRef, U, Trans, useMass);
      }
      RmsData.AddData(nframe, &R, nref);
      // DEBUG
      //mprinterr("%12i %12i %12.4lf\n",nref,nframe,R);
      current++;
    } // END loop over target frames
  } // END loop over reference frames
  progress->Update(max);
  delete progress;

  if (SelectedRef!=NULL) delete SelectedRef;
  if (SelectedTgt!=NULL) delete SelectedTgt;
  return;
}



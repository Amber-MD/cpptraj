#include <cstdio> //sprintf
#include "Action_Rms2d.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"

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
  if (RefTraj!=NULL) delete RefTraj; 
}

// Rms2d::SeparateInit()
/** For use when not part of the action list, i.e. using the Rms2d action
  * to calculate a distance matrix for e.g. clustering.
  */
int Rms2d::SeparateInit(bool nofitIn, char *maskIn) {
  nofit = nofitIn;
  FrameMask.SetMaskString(maskIn);
  return 0;
};

// Rms2d::init()
/** Expected call: rms2d <mask> <refmask> rmsout <filename> [nofit] 
  *                [reftraj <traj> [parm <parmname> | parmindex <#>]] 
  */
// Dataset name will be the last arg checked for. Check order is:
//    1) Keywords
//    2) Masks
//    3) Dataset name
int Rms2d::init() {
  char *mask0, *maskRef, *reftraj;

  // Get keywords
  nofit = actionArgs.hasKey("nofit");
  //useMass = actionArgs.hasKey("mass"); Since parm info no longer stored, useMass is redundant
  rmsdFile = actionArgs.getKeyString("rmsout",NULL);
  reftraj = actionArgs.getKeyString("reftraj",NULL);
  if (reftraj!=NULL) {
    RefParm = PFL->GetParm(actionArgs);
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
  mask0 = actionArgs.getNextMask();
  FrameMask.SetMaskString(mask0);

  // Check if reference will be a series of frames from a trajectory
  if (reftraj!=NULL) {
    // Get RMS mask string for reference trajectory
    maskRef = actionArgs.getNextMask();
    // If no reference mask specified, make same as RMS mask
    if (maskRef==NULL) maskRef=mask0;
    RefMask.SetMaskString(maskRef);
    // Attempt to set up reference trajectory
    RefTraj = new TrajectoryFile();
    if (RefTraj->SetupRead(reftraj, NULL, RefParm)) {
      mprinterr("Error: Rms2d: Could not set up reftraj %s.\n",reftraj);
      return 1;
    }
  }

  mprintf("    RMS2D: Mask [%s]",FrameMask.MaskString());
  if (reftraj!=NULL) {
    // Set up reference trajectory and open
    mprintf(", ref traj %s (mask [%s]) %i frames",RefTraj->TrajName(),
            RefMask.MaskString(),RefTraj->Total_Read_Frames());
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

// Rms2d::setup()
/** Set up frame mask so that only selected atoms in frames will be stored.
  */
int Rms2d::setup() {
  if (FrameMask.SetupMask(currentParm, activeReference, debug)) {
    mprinterr("Error: Rms2d::setup: Could not set up mask [%s] for parm %s\n",
              FrameMask.MaskString(), currentParm->parmName);
    return 1;
  }
  if (FrameMask.None()) {
    mprinterr("Error: Rms2d::setup: No atoms selected for mask [%s], parm %s\n",
              FrameMask.MaskString(), currentParm->parmName);
    return 1;
  }
  return 0;  
}

// Rms2d::action()
/** Store current frame coords according to mask.
  */
int Rms2d::action() {
 
  if (ReferenceCoords.AddCoordsByMask(currentFrame->X, &FrameMask)) return 1;

  return 0;
} 

// Rms2d::Calc2drms()
/** Calculate the RMSD of each frame in ReferenceCoords to each other frame.
  * Since this results in a symmetric matrix use TriangleMatrix to store
  * results.
  */
void Rms2d::Calc2drms(TriangleMatrix *Distances) {
  Frame RefFrame;
  Frame TgtFrame;
  Frame SelectedRef;
  Frame SelectedTgt;
  float *coord;
  double U[9], Trans[6];
  float R;
  int natom_ref, natom_tgt;
  int totalref, max;
 
  totalref = ReferenceCoords.Ncoords();
  Distances->Setup( totalref );
  max = Distances->Nelements();
  mprintf("  RMS2D: Calculating RMSDs between each frame (%i total).\n  ",max);

  // Set up progress Bar
  ProgressBar *progress = new ProgressBar(max);

  // LOOP OVER REFERENCE FRAMES
  int current = 0;
  for (int nref=0; nref < totalref - 1; nref++) {
    progress->Update(current);
    // Get the current reference frame
    coord = ReferenceCoords.Coord(nref, &natom_ref);
    RefFrame.SetupFrameFromCoords( coord, natom_ref );
  
    // LOOP OVER TARGET FRAMES
    for (int nframe=nref+1; nframe < totalref; nframe++) {
      // Get the current target frame
      coord = ReferenceCoords.Coord(nframe, &natom_tgt);
      TgtFrame.SetupFrameFromCoords( coord, natom_tgt );
    
      // Ensure # ref atoms == # tgt atoms
      if (natom_ref != natom_tgt) {
        mprintf("\tWarning: Rms2d: # atoms in ref %i (%i) != # atoms in tgt %i (%i)\n",
                nref+1,natom_ref,nframe+1,natom_tgt);
        mprintf("\t         Assigning this pair RMSD of -1.0\n");
        R = -1.0;
        Distances->AddElement( R );
        continue;
      }
      // Set selected reference atoms - always done since RMS fit modifies SelectedRef
      SelectedRef = RefFrame;
      // Set selected target atoms
      SelectedTgt = TgtFrame;

      // Perform RMS calculation
      if (nofit) {
        R = (float) SelectedTgt.RMSD(&SelectedRef, false);
      } else {
        R = (float) SelectedTgt.RMSD(&SelectedRef, U, Trans, false);
      }
      Distances->AddElement( R );
      // DEBUG
      //mprinterr("%12i %12i %12.4lf\n",nref,nframe,R);
      current++;
    } // END loop over target frames
  } // END loop over reference frames
  progress->Update(max);
  delete progress;
}

// Rms2d::CalcRmsToTraj()
/** Calc RMSD of every frame in reference traj to every frame in 
  * ReferenceCoords.
  */
void Rms2d::CalcRmsToTraj() {
  Frame RefFrame;
  Frame TgtFrame;
  Frame SelectedRef;
  Frame SelectedTgt;
  DataSet *rmsdata;
  float *coord;
  char setname[256];
  double U[9], Trans[6];
  float R;
  int natom_tgt;
  // Set up reference mask for reference parm
  if (RefMask.SetupMask(RefParm, activeReference, debug)) {
    mprinterr("Error: Could not set up reference mask [%s] for parm %s\n",
              RefMask.MaskString(), RefParm->parmName);
    return;
  }
  int natom_ref = RefMask.Nselected;
  // Setup frame for selected reference atoms
  SelectedRef.SetupFrameFromMask(&RefMask, RefParm->mass); 
  RefFrame.SetupFrame(RefParm->natom,RefParm->mass);
  int totalref = RefTraj->Total_Read_Frames();
  int totaltgt = ReferenceCoords.Ncoords();
  int max = totalref * totaltgt;
  mprintf("  RMS2D: Calculating RMSDs between each input frame and each reference\n"); 
  mprintf("         trajectory %s frame (%i total).\n  ",
          RefTraj->TrajName(), max);
  if (RefTraj->BeginTraj(false)) {
    mprinterr("Error: Rms2d: Could not open reference trajectory.\n");
    return;
  }
  // Set up progress Bar
  ProgressBar *progress = new ProgressBar(max);
  // LOOP OVER REFERENCE FRAMES
  int current=0;
  for (int nref=0; nref < totalref; nref++) {
    progress->Update(current);
    // Get the current reference frame from trajectory
    RefTraj->GetNextFrame(RefFrame);
  
    // Set up dataset for this reference frame
    sprintf(setname,"Frame_%i",nref+1);
    rmsdata = RmsData.Add(FLOAT, setname, "Rms2d");
    DFL->Add(rmsdFile,rmsdata);

    // LOOP OVER TARGET FRAMES
    for (int nframe=0; nframe < totaltgt; nframe++) {
      // Get the current target frame
      coord = ReferenceCoords.Coord(nframe, &natom_tgt);
      TgtFrame.SetupFrameFromCoords( coord, natom_tgt );

      // Ensure # ref atoms == # tgt atoms
      if (natom_ref != natom_tgt) {
        mprintf("\tWarning: Rms2d: # atoms in ref %i (%i) != # atoms in tgt %i (%i)\n",
                nref+1,natom_ref,nframe+1,natom_tgt);
        mprintf("\t         Assigning this pair RMSD of -1.0\n");
        R = -1.0;
        RmsData.AddData(nframe, &R, nref);
        continue;
      }

      // Set selected reference atoms - always done since RMS fit modifies SelectedRef
      SelectedRef.SetFrameCoordsFromMask(RefFrame.X, &RefMask);
      // Set selected target atoms
      SelectedTgt = TgtFrame;

      // Perform RMS calculation
      if (nofit) {
        R = (float) SelectedTgt.RMSD(&SelectedRef, false);
      } else {
        R = (float) SelectedTgt.RMSD(&SelectedRef, U, Trans, false);
      }
      RmsData.AddData(nframe, &R, nref);
      // DEBUG
      //mprinterr("%12i %12i %12.4lf\n",nref,nframe,R);
      current++;
    } // END loop over target frames
  } // END loop over reference frames
  progress->Update(max);
  delete progress;
  RefTraj->EndTraj();
}

// Rms2d::print()
/** Perform the rms calculation of each frame to each other frame.
  */
void Rms2d::print() {
  TriangleMatrix *Distances;
  char setname[256];

  if (RefTraj==NULL) {
    Distances = new TriangleMatrix(); 
    Calc2drms(Distances);
    // Convert TriangleMatrix to a dataset.
    // NOTE: This currently uses way more memory than it needs to.
    // TriangleMatrix should just be a dataset that can be output.
    for (int nref=0; nref < ReferenceCoords.Ncoords(); nref++) {
      // Set up dataset for this reference frame
      sprintf(setname,"Frame_%i",nref+1);
      DataSet *rmsdata = RmsData.Add(FLOAT, setname, "Rms2d");
      DFL->Add(rmsdFile,rmsdata);
      for (int nframe=0; nframe < ReferenceCoords.Ncoords(); nframe++) {
        float R = Distances->GetElementF(nref, nframe);
        RmsData.AddData(nframe, &R, nref);
      }
    }
    delete Distances;
  } else
    CalcRmsToTraj();
}


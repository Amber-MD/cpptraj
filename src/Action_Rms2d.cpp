// Rms2d 
#include "Action_Rms2d.h"
#include "CpptrajStdio.h"
#include <cstdio> //sprintf

// CONSTRUCTOR
Rms2d::Rms2d() {
  //fprintf(stderr,"Rms2d Con\n");
  nofit=false;
  useMass=false;
  rmsdFile=NULL;
} 

// DESTRUCTOR
Rms2d::~Rms2d() { }

/*
 * Rms2d::init()
 * Expected call: rms2d <mask> <refmask> [rmsout filename] [mass] [nofit]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Rms2d::init() {
  char *mask0, *maskRef;

  // Get keywords
  nofit = A->hasKey("nofit");
  useMass = A->hasKey("mass");
  rmsdFile = A->getKeyString("rmsout",NULL);
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

  mprintf("    RMS2D: (%s) to (%s)",FrameMask.maskString,RefMask.maskString);
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
  AmberParm *RefParm;
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
  int currentPercent=0;
  int targetPercent=0;
  DataSet *rmsdata;
  char setname[256];

  max = ReferenceFrames.NumFrames() * ReferenceFrames.NumFrames();
  mprintf("  RMS2D: Calculating RMSDs between each frame (%i total).\n  ",max);

  for (int nref=0; nref < ReferenceFrames.NumFrames(); nref++) {
    currentPercent = (current * 100) / max;
    if (currentPercent >= targetPercent) {
      targetPercent+=10;
      mprintf("%2i%% ",currentPercent);
      mflush();
    }
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
    RefFrame = ReferenceFrames.GetFrame( nref );
    // Set up dataset for this reference frame
    sprintf(setname,"Frame_%i",nref);
    rmsdata = RmsData.Add(DOUBLE, setname, "Rms2d");
    DFL->Add(rmsdFile,rmsdata);

    for (int nframe=0; nframe < ReferenceFrames.NumFrames(); nframe++) {
      //progressBar(current, max);
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

      if (nofit) {
        R = SelectedTgt->RMSD(SelectedRef, useMass);
      } else {
        R = SelectedTgt->RMSD(SelectedRef, U, Trans, useMass);
        //F->Translate(Trans);
        //F->Rotate(U);
        //F->Translate(Trans+3);
      }
      RmsData.AddData(nframe, &R, nref);
      // DEBUG
      //mprinterr("%12i %12i %12.4lf\n",nref,nframe,R);
      current++;
    } // END loop over target frames
  } // END loop over reference frames
  //progressBar(current, max);
  mprintf("100%\n");

  if (SelectedRef!=NULL) delete SelectedRef;
  if (SelectedTgt!=NULL) delete SelectedTgt;
  return;
}

/*
 * Rms2d::progressBar()
 * Print progress bar during 2d rms calc
 */
void Rms2d::progressBar(int current, int max) {
  char buffer[128];
  int i,j,percent,numChars;

  // Fraction complete, max 100. 
  percent = (current * 100) / max;

  buffer[0]='{';
  i=1;
  // Set number of characters
  numChars = percent / 2;
  for (j=0; j<numChars; j++) buffer[i++]='+';
  // Fill the rest
  for (; j<50; j++)
    buffer[i++]=' ';
  // Add last character and reset/newline character
  if (current==max) buffer[i-1]='+';
  buffer[i++]='}';
  buffer[i++]='\r';
  if (current==max) {
    buffer[i-1]=' ';
    buffer[i++]='C';
    buffer[i++]='o';
    buffer[i++]='m';
    buffer[i++]='p';
    buffer[i++]='l';
    buffer[i++]='e';
    buffer[i++]='t';
    buffer[i++]='e';
    buffer[i++]='.';
    buffer[i++]='\n';
  }
  // Finish off and print
  buffer[i]='\0';
  mprintf("  %s",buffer);
  // On first frame flush so that progress bar appears
  if (current==0) mflush();
}

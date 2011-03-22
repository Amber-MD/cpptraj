// Rms2d 
#include "Action_Rms2d.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Rms2d::Rms2d() {
  //fprintf(stderr,"Rms2d Con\n");
  nofit=false;
  useMass=false;
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

  ReferenceFrames.Add(F,P);
  
  return 0;
} 



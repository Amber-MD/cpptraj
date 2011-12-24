// Distance
#include <cmath>
#include "Action_Distance.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Distance::Distance() {
  //fprintf(stderr,"Distance Con\n");
  dist=NULL;
  noimage=false;
  useMass=true;
  imageType=0;
} 

// DESTRUCTOR
Distance::~Distance() {
  //fprintf(stderr,"Distance Destructor.\n");
}

// Distance::init()
/// Expected call: distance <name> <mask1> <mask2> [out filename] [geom] [noimage]
// Dataset name will be the last arg checked for. Check order is:
//    1) Keywords
//    2) Masks
//    3) Dataset name
int Distance::init( ) {
  char *mask1, *mask2;
  char *distanceFile;

  // Get Keywords
  noimage = actionArgs.hasKey("noimage");
  useMass = !(actionArgs.hasKey("geom"));
  distanceFile = actionArgs.getKeyString("out",NULL);

  // Get Masks
  mask1 = actionArgs.getNextMask();
  mask2 = actionArgs.getNextMask();
  //fprintf(stdout,"    Mask 1: %s\n",mask1);
  //fprintf(stdout,"    Mask 2: %s\n",mask2);
  if (mask1==NULL || mask2==NULL) {
    mprinterr("Error: Distance::init: Requires 2 masks\n");
    return 1;
  }
  Mask1.SetMaskString(mask1);
  Mask2.SetMaskString(mask2);

  // Dataset to store distances
  dist = DSL->Add(DOUBLE, actionArgs.getNextString(),"Dis");
  if (dist==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(distanceFile,dist);

  mprintf("    DISTANCE: %s to %s",Mask1.MaskString(), Mask2.MaskString());
  if (noimage) 
    mprintf(", non-imaged");
  if (useMass) 
    mprintf(", center of mass");
  else
    mprintf(", geometric center");
  mprintf(".\n");

  return 0;
}

// Distance::setup()
/** Determine what atoms each mask pertains to for the current parm file.
  * Also determine whether imaging should be performed.
  */
int Distance::setup() {

  if (currentParm->SetupIntegerMask( Mask1, activeReference)) return 1;
  if (currentParm->SetupIntegerMask( Mask2, activeReference)) return 1;

  // Check imaging - check box based on prmtop box
  imageType = 0;
  if (!noimage) {
    imageType = (int)currentParm->boxType;
    if (currentParm->boxType==NOBOX && debug>0) {
      mprintf("    Warning: No box info in %s, disabling imaging.\n",currentParm->parmName);
    }
  }

  // Print imaging info for this parm
  mprintf("\t%s (%i atoms) to %s (%i atoms)",Mask1.MaskString(), Mask1.Nselected,
          Mask2.MaskString(),Mask2.Nselected);
  if (imageType > 0)
    mprintf(", imaged");
  else
    mprintf(", imaging off");
  mprintf(".\n");

  if (Mask1.None() || Mask2.None()) {
    mprintf("    Warning: Distance::setup: One or both masks have no atoms.\n");
    return 1;
  }
       
  return 0;  
}

// Distance::action()
int Distance::action() {
  double D, ucell[9], recip[9];

  if (imageType>0) currentFrame->BoxToRecip(ucell,recip);
  D = currentFrame->DIST2(&Mask1, &Mask2, useMass, imageType, ucell, recip);
  D = sqrt(D);

  dist->Add(frameNum, &D);

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,D);
  
  return 0;
} 


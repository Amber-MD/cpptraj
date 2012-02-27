// Angle 
#include "Action_Angle.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG

// CONSTRUCTOR
Angle::Angle() {
  //fprintf(stderr,"Angle Con\n");
  ang=NULL;
  useMass=false;
} 

// Angle::init()
/** Expected call: angle <name> <mask1> <mask2> <mask3> [out filename] [mass]
  */
// Dataset name will be the last arg checked for. Check order is:
//    1) Keywords
//    2) Masks
//    3) Dataset name
int Angle::init() {
  char *mask1, *mask2, *mask3;
  char *angleFile;

  // Get keywords
  angleFile = actionArgs.getKeyString("out",NULL);
  useMass = actionArgs.hasKey("mass");

  // Get Masks
  mask1 = actionArgs.getNextMask();
  mask2 = actionArgs.getNextMask();
  mask3 = actionArgs.getNextMask();
  if (mask1==NULL || mask2==NULL || mask3==NULL) {
    mprinterr("Error: Angle::init: Requires 3 masks\n");
    return 1;
  }
  Mask1.SetMaskString(mask1);
  Mask2.SetMaskString(mask2);
  Mask3.SetMaskString(mask3);

  // Dataset to store angles
  ang = DSL->Add(DOUBLE, actionArgs.getNextString(),"Ang");
  if (ang==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(angleFile,ang);

  mprintf("    ANGLE: [%s]-[%s]-[%s]\n",Mask1.MaskString(),Mask2.MaskString(),Mask3.MaskString());
  if (useMass)
    mprintf("              Using center of mass of atoms in masks.\n");

  return 0;
}

// Angle::setup()
/** Set angle up for this parmtop. Get masks etc.
  */
// currentParm is set in Action::Setup
int Angle::setup() {

  if (currentParm->SetupIntegerMask(Mask1, activeReference)) return 1;
  if (currentParm->SetupIntegerMask(Mask2, activeReference)) return 1;
  if (currentParm->SetupIntegerMask(Mask3, activeReference)) return 1;
  mprintf("\t%s (%i atoms)\n",Mask1.MaskString(),Mask1.Nselected);
  mprintf("\t%s (%i atoms)\n",Mask2.MaskString(),Mask2.Nselected);
  mprintf("\t%s (%i atoms)\n",Mask3.MaskString(),Mask3.Nselected);
  if (Mask1.None() || Mask2.None() || Mask3.None()) {
    mprintf("    Error: Angle::setup: One or more masks contain 0 atoms.\n");
    return 1;
  }


  return 0;  
}

// Angle::action()
int Angle::action() {
  double Ang;

  Ang=currentFrame->ANGLE(&Mask1,&Mask2,&Mask3,useMass);

  Ang *= RADDEG;

  ang->Add(frameNum, &Ang);

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,Ang);
  
  return 0;
} 


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

// DESTRUCTOR
Angle::~Angle() { }

/* Angle::init()
 * Expected call: angle <name> <mask1> <mask2> <mask3> [out filename] [mass]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Angle::init() {
  char *mask1, *mask2, *mask3;
  char *angleFile;

  // Get keywords
  angleFile = A->getKeyString("out",NULL);
  useMass = A->hasKey("mass");

  // Get Masks
  mask1 = A->getNextMask();
  mask2 = A->getNextMask();
  mask3 = A->getNextMask();
  if (mask1==NULL || mask2==NULL || mask3==NULL) {
    mprintf("    Error: Angle::init: Requires 3 masks\n");
    return 1;
  }
  Mask1.SetMaskString(mask1);
  Mask2.SetMaskString(mask2);
  Mask3.SetMaskString(mask3);

  // Dataset to store angles
  ang = DSL->Add(DOUBLE, A->getNextString(),"Ang");
  if (ang==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(angleFile,ang);

  mprintf("    ANGLE: %s-%s-%s\n",Mask1.maskString,Mask2.maskString,Mask3.maskString);
  if (useMass)
    mprintf("              Using center of mass of atoms in masks.\n");

  return 0;
}

/* Angle::setup()
 * Set angle up for this parmtop. Get masks etc.
 * P is set in Action::Setup
 */
int Angle::setup() {

  if ( Mask1.SetupMask(P,debug) ) return 1;
  if ( Mask2.SetupMask(P,debug) ) return 1;
  if ( Mask3.SetupMask(P,debug) ) return 1;
  if (Mask1.None() || Mask2.None() || Mask3.None()) {
    mprintf("    Error: Angle::setup: One or more masks contain 0 atoms.\n");
    return 1;
  }

  return 0;  
}

/* Angle::action()
 */
int Angle::action() {
  double Ang;

  Ang=F->ANGLE(&Mask1,&Mask2,&Mask3,useMass);

  Ang *= RADDEG;

  ang->Add(currentFrame, &Ang);

  //fprintf(outfile,"%10i %10.4lf\n",currentFrame,Ang);
  
  return 0;
} 



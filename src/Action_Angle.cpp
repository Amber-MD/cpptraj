// Angle 
#include "Action_Angle.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG

// CONSTRUCTOR
Angle::Angle() :
  ang_(NULL)
{
  //fprintf(stderr,"Angle Con\n");
  useMass_ = false;
} 

// Angle::init()
/** Expected call: angle <name> <mask1> <mask2> <mask3> [out filename] [mass]
  */
// Dataset name will be the last arg checked for. Check order is:
//    1) Keywords
//    2) Masks
//    3) Dataset name
int Angle::init() {
  // Get keywords
  char *angleFile = actionArgs.getKeyString("out",NULL);
  useMass_ = actionArgs.hasKey("mass");

  // Get Masks
  char *mask1 = actionArgs.getNextMask();
  char *mask2 = actionArgs.getNextMask();
  char *mask3 = actionArgs.getNextMask();
  if (mask1==NULL || mask2==NULL || mask3==NULL) {
    mprinterr("Error: Angle::init: Requires 3 masks\n");
    return 1;
  }
  Mask1_.SetMaskString(mask1);
  Mask2_.SetMaskString(mask2);
  Mask3_.SetMaskString(mask3);

  // Dataset to store angles
  ang_ = DSL->Add(DataSet::DOUBLE, actionArgs.getNextString(),"Ang");
  if (ang_==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(angleFile, ang_);

  mprintf("    ANGLE: [%s]-[%s]-[%s]\n",Mask1_.MaskString(), Mask2_.MaskString(), 
          Mask3_.MaskString());
  if (useMass_)
    mprintf("              Using center of mass of atoms in masks.\n");

  return 0;
}

// Angle::setup()
/** Set angle up for this parmtop. Get masks etc.
  */
// currentParm is set in Action::Setup
int Angle::setup() {

  if (currentParm->SetupIntegerMask(Mask1_)) return 1;
  if (currentParm->SetupIntegerMask(Mask2_)) return 1;
  if (currentParm->SetupIntegerMask(Mask3_)) return 1;
  mprintf("\t%s (%i atoms)\n",Mask1_.MaskString(), Mask1_.Nselected());
  mprintf("\t%s (%i atoms)\n",Mask2_.MaskString(), Mask2_.Nselected());
  mprintf("\t%s (%i atoms)\n",Mask3_.MaskString(), Mask3_.Nselected());
  if (Mask1_.None() || Mask2_.None() || Mask3_.None()) {
    mprintf("    Error: Angle::setup: One or more masks contain 0 atoms.\n");
    return 1;
  }


  return 0;  
}

// Angle::action()
int Angle::action() {
  double aval = currentFrame->ANGLE(Mask1_, Mask2_, Mask3_, useMass_);

  aval *= RADDEG;

  ang_->Add(frameNum, &aval);

  return 0;
}


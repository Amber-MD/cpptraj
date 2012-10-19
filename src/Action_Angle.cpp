// Action_Angle 
#include "Action_Angle.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG
#include "TorsionRoutines.h"

// CONSTRUCTOR
Action_Angle::Action_Angle() :
  ang_(NULL),
  useMass_(false)
{ } 

void Action_Angle::Help() {

}

// Action_Angle::init()
/** Expected call: angle <name> <mask1> <mask2> <mask3> [out filename] [mass]
  */
// Dataset name will be the last arg checked for. Check order is:
//    1) Keywords
//    2) Masks
//    3) Dataset name
int Action_Angle::init() {
  // Get keywords
  ArgList::ConstArg angleFile = actionArgs.getKeyString("out");
  useMass_ = actionArgs.hasKey("mass");

  // Get Masks
  ArgList::ConstArg mask1 = actionArgs.getNextMask();
  ArgList::ConstArg mask2 = actionArgs.getNextMask();
  ArgList::ConstArg mask3 = actionArgs.getNextMask();
  if (mask1==NULL || mask2==NULL || mask3==NULL) {
    mprinterr("Error: angle: Requires 3 masks\n");
    return 1;
  }
  Mask1_.SetMaskString(mask1);
  Mask2_.SetMaskString(mask2);
  Mask3_.SetMaskString(mask3);

  // Dataset to store angles
  ang_ = DSL->Add(DataSet::DOUBLE, actionArgs.getNextString(),"Ang");
  if (ang_==NULL) return 1;
  ang_->SetScalar( DataSet::M_ANGLE );
  // Add dataset to data file list
  DFL->AddSetToFile(angleFile, ang_);

  mprintf("    ANGLE: [%s]-[%s]-[%s]\n",Mask1_.MaskString(), Mask2_.MaskString(), 
          Mask3_.MaskString());
  if (useMass_)
    mprintf("              Using center of mass of atoms in masks.\n");

  return 0;
}

// Action_Angle::setup()
/** Set angle up for this parmtop. Get masks etc.
  */
// currentParm is set in Action::Setup
int Action_Angle::setup() {

  if (currentParm->SetupIntegerMask(Mask1_)) return 1;
  if (currentParm->SetupIntegerMask(Mask2_)) return 1;
  if (currentParm->SetupIntegerMask(Mask3_)) return 1;
  Mask1_.MaskInfo();
  Mask2_.MaskInfo();
  Mask3_.MaskInfo();
  if (Mask1_.None() || Mask2_.None() || Mask3_.None()) {
    mprintf("Warning: angle: One or more masks contain 0 atoms.\n");
    return 1;
  }

  return 0;  
}

// Action_Angle::action()
int Action_Angle::action() {
  Vec3 a1, a2, a3;
  if (useMass_) {
    a1 = currentFrame->VCenterOfMass( Mask1_ );
    a2 = currentFrame->VCenterOfMass( Mask2_ );
    a3 = currentFrame->VCenterOfMass( Mask3_ );
  } else {
    a1 = currentFrame->VGeometricCenter( Mask1_ );
    a2 = currentFrame->VGeometricCenter( Mask2_ );
    a3 = currentFrame->VGeometricCenter( Mask3_ );
  }
  double aval = CalcAngle( a1.Dptr(), a2.Dptr(), a3.Dptr() );

  aval *= RADDEG;

  ang_->Add(frameNum, &aval);

  return 0;
}


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
  mprintf("angle [<name>] <mask1> <mask2> <mask3> [out filename] [mass]\n");
  mprintf("\tCalculate the angle between atoms in masks 1-3\n");
}

// Action_Angle::init()
Action::RetType Action_Angle::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  ArgList::ConstArg angleFile = actionArgs.getKeyString("out");
  useMass_ = actionArgs.hasKey("mass");

  // Get Masks
  ArgList::ConstArg mask1 = actionArgs.getNextMask();
  ArgList::ConstArg mask2 = actionArgs.getNextMask();
  ArgList::ConstArg mask3 = actionArgs.getNextMask();
  if (mask1==NULL || mask2==NULL || mask3==NULL) {
    mprinterr("Error: angle: Requires 3 masks\n");
    return Action::ERR;
  }
  Mask1_.SetMaskString(mask1);
  Mask2_.SetMaskString(mask2);
  Mask3_.SetMaskString(mask3);

  // Dataset to store angles
  ang_ = DSL->Add(DataSet::DOUBLE, actionArgs.getNextString(),"Ang");
  if (ang_==NULL) return Action::ERR;
  ang_->SetScalar( DataSet::M_ANGLE );
  // Add dataset to data file list
  DFL->AddSetToFile(angleFile, ang_);

  mprintf("    ANGLE: [%s]-[%s]-[%s]\n",Mask1_.MaskString(), Mask2_.MaskString(), 
          Mask3_.MaskString());
  if (useMass_)
    mprintf("              Using center of mass of atoms in masks.\n");

  return Action::OK;
}

// Action_Angle::setup()
/** Set angle up for this parmtop. Get masks etc.
  */
// currentParm is set in Action::Setup
Action::RetType Action_Angle::Setup(Topology* currentParm, Topology** parmAddress) {

  if (currentParm->SetupIntegerMask(Mask1_)) return Action::ERR;
  if (currentParm->SetupIntegerMask(Mask2_)) return Action::ERR;
  if (currentParm->SetupIntegerMask(Mask3_)) return Action::ERR;
  Mask1_.MaskInfo();
  Mask2_.MaskInfo();
  Mask3_.MaskInfo();
  if (Mask1_.None() || Mask2_.None() || Mask3_.None()) {
    mprintf("Warning: angle: One or more masks contain 0 atoms.\n");
    return Action::ERR;
  }

  return Action::OK;  
}

// Action_Angle::action()
Action::RetType Action_Angle::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
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

  return Action::OK;
}


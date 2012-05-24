#include "Action_Scale.h"
#include "CpptrajStdio.h"

Action_Scale::Action_Scale() :
  sx_(1),
  sy_(1),
  sz_(1)
{}

int Action_Scale::init() {
  sx_ = actionArgs.getKeyDouble("x", 1);
  sy_ = actionArgs.getKeyDouble("y", 1);
  sz_ = actionArgs.getKeyDouble("z", 1);
  mask_.SetMaskString( actionArgs.getNextMask() );

  mprintf("    SCALE coordinates: X by %.3f, Y by %.3f, Z by %.3f\n", sx_, sy_, sz_);
  mprintf("                       Mask is [%s]\n", mask_.MaskString());

  return 0;
}

int Action_Scale::setup() {
  if ( currentParm->SetupIntegerMask( mask_ ) ) return 1;
  if ( mask_.None() ) {
    mprintf("Warning: scale: No atoms selected.\n");
    return 1;
  }
  return 0;
}

int Action_Scale::action() {
  
  currentFrame->SCALE(mask_, sx_, sy_, sz_);

  return 0;
}


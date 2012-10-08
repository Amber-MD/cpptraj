#include "Action_Rotate.h"
#include "CpptrajStdio.h"
#include "vectormath.h"
#include "Constants.h"

// CONSTRUCTOR
Action_Rotate::Action_Rotate() {
  RotMatrix_[0] = 0.0;
  RotMatrix_[1] = 0.0;
  RotMatrix_[2] = 0.0;
  RotMatrix_[3] = 0.0;
  RotMatrix_[4] = 0.0;
  RotMatrix_[5] = 0.0;
  RotMatrix_[6] = 0.0;
  RotMatrix_[7] = 0.0;
  RotMatrix_[8] = 0.0;
}

/** Usage: rotate [<mask>] [x <xdeg>] [y <ydeg>] [z <zdeg>]
  */
int Action_Rotate::init() {
  double xrot = actionArgs.getKeyDouble("x",0.0);
  double yrot = actionArgs.getKeyDouble("y",0.0);
  double zrot = actionArgs.getKeyDouble("z",0.0);
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  // Calc rotation matrix
  calcRotationMatrix( RotMatrix_, xrot * DEGRAD, yrot * DEGRAD, zrot * DEGRAD );

  mprintf("    ROTATE: Rotating atoms in mask %s\n", mask_.MaskString());
  mprintf("\t%f degrees around X, %f degrees around Y, %f degrees around Z\n",
          xrot, yrot, zrot);
  return 0;
};

int Action_Rotate::setup() {
  if ( currentParm->SetupIntegerMask( mask_ ) ) return 1;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprinterr("Error: rotate: No atoms selected.\n");
    return 1;
  }
  return 0;
}

int Action_Rotate::action() {
  currentFrame->Rotate(RotMatrix_, mask_);
  return 0;
}


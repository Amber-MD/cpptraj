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

void Action_Rotate::Help() {

}

/** Usage: rotate [<mask>] [x <xdeg>] [y <ydeg>] [z <zdeg>]
  */
Action::RetType Action_Rotate::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  double xrot = actionArgs.getKeyDouble("x",0.0);
  double yrot = actionArgs.getKeyDouble("y",0.0);
  double zrot = actionArgs.getKeyDouble("z",0.0);
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  // Calc rotation matrix
  calcRotationMatrix( RotMatrix_, xrot * DEGRAD, yrot * DEGRAD, zrot * DEGRAD );

  mprintf("    ROTATE: Rotating atoms in mask %s\n", mask_.MaskString());
  mprintf("\t%f degrees around X, %f degrees around Y, %f degrees around Z\n",
          xrot, yrot, zrot);
  return Action::OK;
};

Action::RetType Action_Rotate::Setup(Topology* currentParm, Topology** parmAddress) {
  if ( currentParm->SetupIntegerMask( mask_ ) ) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprinterr("Error: rotate: No atoms selected.\n");
    return Action::ERR;
  }
  return Action::OK;
}

Action::RetType Action_Rotate::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  currentFrame->Rotate(RotMatrix_, mask_);
  return Action::OK;
}


#include "Action_Rotate.h"
#include "CpptrajStdio.h"
#include "Constants.h"

// CONSTRUCTOR
Action_Rotate::Action_Rotate() :
  rmatrices_(0), delta_(0.0), mode_(ROTATE), inverse_(false) { }

void Action_Rotate::Help() const {
  mprintf("\t[<mask>] { [x <xdeg>] [y <ydeg>] [z <zdeg>]  |\n"
          "\t           axis0 <mask0> axis1 <mask1> <deg> |\n"
          "\t           usedata <set name> [inverse] }\n"
          "  Rotate atoms in <mask> either around the x, y, and/or z axes, around the\n"
          "  the axis defined by <mask0> to <mask1>, or using rotation matrices in\n"
          "  the specified data set.\n");
}

Action::RetType Action_Rotate::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  double xrot = 0.0, yrot = 0.0, zrot = 0.0;
  std::string dsname = actionArgs.GetStringKey("usedata");
  std::string axis = actionArgs.GetStringKey("axis0");
  if (!dsname.empty()) {
    inverse_ = actionArgs.hasKey("inverse");
    // Check if DataSet exists
    rmatrices_ = (DataSet_Mat3x3*)init.DSL().FindSetOfType( dsname, DataSet::MAT3X3 );
    if (rmatrices_ == 0) {
      mprinterr("Error: No 3x3 matrices data set '%s'\n", dsname.c_str());
      return Action::ERR;
    }
    mode_ = DATASET;
  } else if (!axis.empty()) {
    // Get axis definition
    if (axis0_.SetMaskString( axis )) return Action::ERR;
    axis = actionArgs.GetStringKey("axis1");
    if (axis.empty()) {
      mprinterr("Error: 'axis1' must be specified if 'axis0' is.\n");
      return Action::ERR;
    }
    if (axis1_.SetMaskString( axis )) return Action::ERR;
    delta_ = actionArgs.getNextDouble(0.0);
    if ( !(delta_ > 0.0) && !(delta_ < 0.0) ) {
      mprinterr("Error: Must specify non-zero rotation.\n");
      return Action::ERR;
    }
    mode_ = AXIS;
  } else {
    // Calc rotation matrix
    xrot = actionArgs.getKeyDouble("x",0.0);
    yrot = actionArgs.getKeyDouble("y",0.0);
    zrot = actionArgs.getKeyDouble("z",0.0);
    RotMatrix_.CalcRotationMatrix( xrot * Constants::DEGRAD, 
                                   yrot * Constants::DEGRAD, 
                                   zrot * Constants::DEGRAD );
  }
  // Get mask
  if (mask_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;

  mprintf("    ROTATE: Rotating atoms in mask %s\n", mask_.MaskString());
  switch (mode_) {
    case ROTATE:
      mprintf("\t%f degrees around X, %f degrees around Y, %f degrees around Z\n",
              xrot, yrot, zrot);
      break;
    case DATASET:
      mprintf("\tUsing rotation matrices from set '%s'\n", rmatrices_->legend());
      if (inverse_) mprintf("\tPerforming inverse rotation.\n");
      break;
    case AXIS:
      mprintf("\t%f degrees around axis defined by '%s' and '%s'\n",
              delta_, axis0_.MaskString(), axis1_.MaskString());
      delta_ *= Constants::DEGRAD;
      break;
  }
  return Action::OK;
}

Action::RetType Action_Rotate::Setup(ActionSetup& setup) {
  if ( setup.Top().SetupIntegerMask( mask_ ) ) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprintf("Warning: No atoms selected.\n");
    return Action::SKIP;
  }
  if (mode_ == AXIS) {
    if ( setup.Top().SetupIntegerMask( axis0_ ) ||
         setup.Top().SetupIntegerMask( axis1_ ) )
      return Action::ERR;
    axis0_.MaskInfo();
    axis1_.MaskInfo();
    if (axis0_.None() || axis1_.None()) {
      mprintf("Warning: Not enough atoms selected to define axis.\n");
      return Action::SKIP;
    }
  }
  Action::CheckImageRotationWarning(setup, "the rotation");
  return Action::OK;
}

Action::RetType Action_Rotate::DoAction(int frameNum, ActionFrame& frm) {
  switch (mode_) {
    case ROTATE : frm.ModifyFrm().Rotate(RotMatrix_, mask_); break;
    case DATASET:
      if (frm.TrajoutNum() >= (int)rmatrices_->Size()) {
        mprintf("Warning: Frame %i out of range for set '%s'\n",
                frm.TrajoutNum()+1, rmatrices_->legend());
        return Action::ERR;
      }
      if (inverse_)
        frm.ModifyFrm().InverseRotate((*rmatrices_)[frm.TrajoutNum()], mask_);
      else
        frm.ModifyFrm().Rotate((*rmatrices_)[frm.TrajoutNum()], mask_);
      break;
    case AXIS   :
      Vec3 a0 = frm.Frm().VCenterOfMass(axis0_);
      Vec3 axisOfRotation = frm.ModifyFrm().SetAxisOfRotation( a0,
                                                               frm.Frm().VCenterOfMass(axis1_) );
      RotMatrix_.CalcRotationMatrix(axisOfRotation, delta_);
      frm.ModifyFrm().Rotate(RotMatrix_, mask_);
      // SetAxisOfRotation moves a0 to center; move back.
      frm.ModifyFrm().Translate( a0 ); 
      break;
  }
  return Action::MODIFY_COORDS;
}

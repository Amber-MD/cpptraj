#include "Action_Rotate.h"
#include "CpptrajStdio.h"
#include "Constants.h"

// CONSTRUCTOR
Action_Rotate::Action_Rotate() : rmatrices_(0), inverse_(false) { }

void Action_Rotate::Help() const {
  mprintf("\t[<mask>] {[x <xdeg>] [y <ydeg>] [z <zdeg>] | usedata <set name>} [inverse]\n"
          "  Rotate atoms in <mask> around x, y, and/or z axes.\n");
}

Action::RetType Action_Rotate::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  double xrot = 0.0, yrot = 0.0, zrot = 0.0;
  inverse_ = actionArgs.hasKey("inverse");
  std::string dsname = actionArgs.GetStringKey("usedata");
  if (dsname.empty()) {
    // Calc rotation matrix
    xrot = actionArgs.getKeyDouble("x",0.0);
    yrot = actionArgs.getKeyDouble("y",0.0);
    zrot = actionArgs.getKeyDouble("z",0.0);
    RotMatrix_.CalcRotationMatrix( xrot * Constants::DEGRAD, 
                                   yrot * Constants::DEGRAD, 
                                   zrot * Constants::DEGRAD );
  } else {
    // Check if DataSet exists
    rmatrices_ = (DataSet_Mat3x3*)init.DSL().FindSetOfType( dsname, DataSet::MAT3X3 );
    if (rmatrices_ == 0) {
      mprinterr("Error: No 3x3 matrices data set '%s'\n", dsname.c_str());
      return Action::ERR;
    }
  }
  // Get mask
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("    ROTATE: Rotating atoms in mask %s\n", mask_.MaskString());
  if (rmatrices_ == 0)
    mprintf("\t%f degrees around X, %f degrees around Y, %f degrees around Z\n",
            xrot, yrot, zrot);
  else
    mprintf("\tUsing rotation matrices from set '%s'\n", rmatrices_->legend());
  if (inverse_) mprintf("\tPerforming inverse rotation.\n");
  return Action::OK;
};

Action::RetType Action_Rotate::Setup(ActionSetup& setup) {
  if ( setup.Top().SetupIntegerMask( mask_ ) ) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprintf("Warning: No atoms selected.\n");
    return Action::SKIP;
  }
  return Action::OK;
}

Action::RetType Action_Rotate::DoAction(int frameNum, ActionFrame& frm) {
  if (rmatrices_ == 0) {
    if (inverse_)
      frm.ModifyFrm().InverseRotate(RotMatrix_, mask_);
    else
      frm.ModifyFrm().Rotate(RotMatrix_, mask_);
  } else {
    if (frameNum >= (int)rmatrices_->Size()) {
      mprintf("Warning: Frame %i out of range for set '%s'\n", frameNum+1, rmatrices_->legend());
      return Action::ERR;
    }
    if (inverse_)
      frm.ModifyFrm().InverseRotate((*rmatrices_)[frameNum], mask_);
    else
      frm.ModifyFrm().Rotate((*rmatrices_)[frameNum], mask_);
  }
  return Action::MODIFY_COORDS;
}

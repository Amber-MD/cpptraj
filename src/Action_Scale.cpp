#include "Action_Scale.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Scale::Action_Scale() :
  sx_(1),
  sy_(1),
  sz_(1)
{}

void Action_Scale::Help() const {
  mprintf("\t[x <sx>] [y <sy>] [z <sz>] [<mask>]\n"
          "\tScale the position of atoms in <mask>\n");
}

// Action_Scale::init()
Action::RetType Action_Scale::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  sx_ = actionArgs.getKeyDouble("x", 1);
  sy_ = actionArgs.getKeyDouble("y", 1);
  sz_ = actionArgs.getKeyDouble("z", 1);
  if (mask_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;

  mprintf("    SCALE coordinates: X by %.3f, Y by %.3f, Z by %.3f\n", sx_, sy_, sz_);
  mprintf("                       Mask is [%s]\n", mask_.MaskString());

  return Action::OK;
}

// Action_Scale::setup()
Action::RetType Action_Scale::Setup(ActionSetup& setup) {
  if ( setup.Top().SetupIntegerMask( mask_ ) ) return Action::ERR;
  if ( mask_.None() ) {
    mprintf("Warning: scale: No atoms selected.\n");
    return Action::SKIP;
  }
  return Action::OK;
}

// Action_Scale::action()
Action::RetType Action_Scale::DoAction(int frameNum, ActionFrame& frm) {
  frm.ModifyFrm().Scale(mask_, sx_, sy_, sz_);
  return Action::MODIFY_COORDS;
}

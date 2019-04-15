#include "Action_Translate.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Translate::Action_Translate() { }

void Action_Translate::Help() const {
  mprintf("\t[<mask>] [x <dx>] [y <dy>] [z <dz>]\n"
          "  Translate atoms in <mask>\n");
}

Action::RetType Action_Translate::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  double x = actionArgs.getKeyDouble("x",0.0);
  double y = actionArgs.getKeyDouble("y",0.0);
  double z = actionArgs.getKeyDouble("z",0.0);
  Trans_.SetVec(x, y, z);
  if (mask_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;

  mprintf("    TRANSLATE: Translating atoms in mask %s\n", mask_.MaskString());
  mprintf("\t%f Ang. in X, %f Ang. in Y, %f Ang. in Z\n",
          Trans_[0], Trans_[1], Trans_[2]);
  return Action::OK;
};

Action::RetType Action_Translate::Setup(ActionSetup& setup) {
  if ( setup.Top().SetupIntegerMask( mask_ ) ) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprintf("Warning: translate: No atoms selected.\n");
    return Action::SKIP;
  }
  return Action::OK;
}

Action::RetType Action_Translate::DoAction(int frameNum, ActionFrame& frm) {
  for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
    frm.ModifyFrm().Translate(Trans_, *atom);
  return Action::MODIFY_COORDS;
}

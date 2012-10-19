#include "Action_Translate.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Translate::Action_Translate() {
  Trans_[0] = 0.0;
  Trans_[1] = 0.0;
  Trans_[2] = 0.0;
}

void Action_Translate::Help() {

}

/** Usage: trans [<mask>] [x <dx>] [y <dy>] [z <dz>]
  */
int Action_Translate::init() {
  Trans_[0] = actionArgs.getKeyDouble("x",0.0);
  Trans_[1] = actionArgs.getKeyDouble("y",0.0);
  Trans_[2] = actionArgs.getKeyDouble("z",0.0);
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("    TRANSLATE: Translating atoms in mask %s\n", mask_.MaskString());
  mprintf("\t%f Ang. in X, %f Ang. in Y, %f Ang. in Z\n",
          Trans_[0], Trans_[1], Trans_[2]);
  return 0;
};

int Action_Translate::setup() {
  if ( currentParm->SetupIntegerMask( mask_ ) ) return 1;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprinterr("Error: translate: No atoms selected.\n");
    return 1;
  }
  return 0;
}

int Action_Translate::action() {
  for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
    currentFrame->Translate(Trans_, *atom);
  return 0;
}


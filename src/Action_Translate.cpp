#include "Action_Translate.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Translate::Action_Translate() :
  Trans_(0.0),
  toPoint_(false),
  useMass_(false)
{ }

void Action_Translate::Help() const {
  mprintf("\t[<mask>] {[x <dx>] [y <dy>] [z <dz>] | topoint <x>,<y>,<z> [mass]}\n"
          "  Translate atoms in <mask> by <dx>, <dy>, and/or <dz> Ang.\n"
          "  If 'topoint' is specified, translate atoms in <mask> to\n"
          "  the specified coordinates (also in Ang.).\n");
}

/** Init. */
Action::RetType Action_Translate::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  std::string toPointArg = actionArgs.GetStringKey("topoint");
  if (toPointArg.empty()) {
    double x = actionArgs.getKeyDouble("x",0.0);
    double y = actionArgs.getKeyDouble("y",0.0);
    double z = actionArgs.getKeyDouble("z",0.0);
    Trans_.SetVec(x, y, z);
  } else {
    ArgList toPointCoords(toPointArg, ",");
    if (toPointCoords.Nargs() != 3) {
      mprinterr("Error: Expected comma-separated list of 3 coordinates for 'topoint'.\n");
      return Action::ERR;
    }
    useMass_ = actionArgs.hasKey("mass");
    toPoint_ = true;
    double x = toPointCoords.getNextDouble(0);
    double y = toPointCoords.getNextDouble(0);
    double z = toPointCoords.getNextDouble(0);
    Trans_.SetVec(x, y, z);
  }
  if (mask_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;

  mprintf("    TRANSLATE: Translating atoms in mask %s\n", mask_.MaskString());
  if (toPoint_) {
    mprintf("\tTarget coords: %f %f %f\n", Trans_[0], Trans_[1], Trans_[2]);
    if (useMass_)
      mprintf("\tWeighting coordinates by mass.\n");
  } else {
    mprintf("\t%f Ang. in X, %f Ang. in Y, %f Ang. in Z\n",
            Trans_[0], Trans_[1], Trans_[2]);
  }
  return Action::OK;
};

/** Setup */
Action::RetType Action_Translate::Setup(ActionSetup& setup) {
  if ( setup.Top().SetupIntegerMask( mask_ ) ) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprintf("Warning: translate: No atoms selected.\n");
    return Action::SKIP;
  }
  return Action::OK;
}

/** DoAction. */
Action::RetType Action_Translate::DoAction(int frameNum, ActionFrame& frm) {
  if (toPoint_) {
    // Calc center of selected coords
    Vec3 ctr;
    if (useMass_)
      ctr = frm.Frm().VCenterOfMass( mask_ );
    else
      ctr = frm.Frm().VGeometricCenter( mask_ );
    // Calc translation to target coord
    Vec3 toTgt = Trans_ - ctr;
    // Translate the atoms
    for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
      frm.ModifyFrm().Translate(toTgt, *atom);
  } else {
    for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
      frm.ModifyFrm().Translate(Trans_, *atom);
  }
  return Action::MODIFY_COORDS;
}

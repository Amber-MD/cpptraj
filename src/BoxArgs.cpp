#include "BoxArgs.h"
#include "ArgList.h"
#include "CpptrajStdio.h"
#include "Box.h"

/** CONSTRUCTOR */
BoxArgs::BoxArgs() {
  for (int i = 0; i < 6; i++) {
    xyzabg_[i] = 0;
    setVar_[i] = false;
  }
}

/** Keywords for SetBoxArgs() */
const char* BoxArgs::Keywords_XyzAbg() {
  return "[x <xval>] [y <yval>] [z <zval>] [alpha <a>] [beta <b>] [gamma <g>]";
}

/** Keyword for truncated oct. */
const char* BoxArgs::Keywords_TruncOct() {
  return "truncoct [x <length>]";
}

/** Set variables from argument list. */
int BoxArgs::SetBoxArgs(ArgList& actionArgs) {
  // TODO check for bad args?
  if (actionArgs.Contains("x")) { xyzabg_[0] = actionArgs.getKeyDouble("x", 0.0); setVar_[0] = true; }
  if (actionArgs.Contains("y")) { xyzabg_[1] = actionArgs.getKeyDouble("y", 0.0); setVar_[1] = true; }
  if (actionArgs.Contains("z")) { xyzabg_[2] = actionArgs.getKeyDouble("z", 0.0); setVar_[2] = true; }
  if (actionArgs.Contains("alpha")) { xyzabg_[3] = actionArgs.getKeyDouble("alpha", 0.0); setVar_[3] = true; }
  if (actionArgs.Contains("beta"))  { xyzabg_[4] = actionArgs.getKeyDouble("beta",  0.0); setVar_[4] = true; }
  if (actionArgs.Contains("gamma")) { xyzabg_[5] = actionArgs.getKeyDouble("gamma", 0.0); setVar_[5] = true; }
  if (actionArgs.hasKey("truncoct")) {
    xyzabg_[3] = Box::TruncatedOctAngle();
    xyzabg_[4] = xyzabg_[3];
    xyzabg_[5] = xyzabg_[3];
    setVar_[3] = true;
    setVar_[4] = true;
    setVar_[5] = true;
    // All lengths need to be the same
    if (setVar_[1]) mprintf("Warning: Only 'x' used for 'truncoct'\n");
    if (setVar_[2]) mprintf("Warning: Only 'x' used for 'truncoct'\n");
    if (setVar_[0]) {
      xyzabg_[1] = xyzabg_[0];
      xyzabg_[2] = xyzabg_[0];
    }
    setVar_[1] = false;
    setVar_[2] = false;
  }
  return 0;
}

/** Set all angles to given value. */
void BoxArgs::SetAngles(double angleIn) {
  for (int i = 3; i < 6; i++)
  {
    xyzabg_[i] = angleIn;
    setVar_[i] = true;
  }
}

/** Set all lengths to given value. */
void BoxArgs::SetLengths(double lengthIn) {
  for (int i = 0; i < 3; i++)
  {
    xyzabg_[i] = lengthIn;
    setVar_[i] = true;
  }
}

/** Print set values of XYZ ABG array to STDOUT. */
void BoxArgs::PrintXyzAbg() const {
  if (setVar_[0]) mprintf(" X=%.3f", xyzabg_[0]);
  if (setVar_[1]) mprintf(" Y=%.3f", xyzabg_[1]);
  if (setVar_[2]) mprintf(" Z=%.3f", xyzabg_[2]);
  if (setVar_[3]) mprintf(" A=%.3f", xyzabg_[3]);
  if (setVar_[4]) mprintf(" B=%.3f", xyzabg_[4]);
  if (setVar_[5]) mprintf(" G=%.3f", xyzabg_[5]);
  mprintf("\n");
}

/** Set any values in xyz abg array not already set with info from box. */
void BoxArgs::SetMissingInfo(Box const& boxIn) {
  for (int i = 0; i < 6; i++) {
    if (!setVar_[i]) xyzabg_[i] = boxIn.Param((Box::ParamType)i);
  }
}

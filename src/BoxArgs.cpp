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
    if (setVar_[1]) { mprintf("Warning: Only 'x' used for 'truncoct'\n"); setVar_[1] = false; }
    if (setVar_[2]) { mprintf("Warning: Only 'x' used for 'truncoct'\n"); setVar_[2] = false; }
    if (setVar_[0]) {
      xyzabg_[1] = xyzabg_[0];
      xyzabg_[2] = xyzabg_[0];
      setVar_[1] = true;
      setVar_[2] = true;
    }
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

/** Used to set an empty param1 from the other params. */
int BoxArgs::SetEmptyInfo(const char* str1, double& param1,
                          const char* str2, double param2,
                          const char* str3, double param3)
{
  mprintf("Warning: Box %s is empty.", str1);
  if (param2 > 0) {
    mprintf(" Setting from %s (%f)\n", str2, param2);
    param1 = param2;
    return 0;
  } else if (param3 > 0) {
    mprintf(" Setting from %s (%f)\n", str3, param3);
    param1 = param3;
    return 0;
  }
  mprintf("\n");
  mprinterr("Error: Nothing available to set box %s\n", str1);
  return 1;
}

/** Set any values in xyz abg array not already set with info from box. */
int BoxArgs::SetMissingInfo(Box const& boxIn) {
  for (int i = 0; i < 6; i++) {
    if (!setVar_[i]) xyzabg_[i] = boxIn.Param((Box::ParamType)i);
    //if (!setVar_[i]) mprintf("DEBUG: SetMissingInfo param %i boxIn= %12.4f\n", i, boxIn.Param((Box::ParamType)i));
  }
  if (xyzabg_[Box::X] <= 0) {
    if (SetEmptyInfo("X", xyzabg_[Box::X], "Y", xyzabg_[Box::Y], "Z", xyzabg_[Box::Z])) return 1;
  }
  if (xyzabg_[Box::Y] <= 0) {
    if (SetEmptyInfo("Y", xyzabg_[Box::Y], "X", xyzabg_[Box::X], "Z", xyzabg_[Box::Z])) return 1;
  }
  if (xyzabg_[Box::Z] <= 0) {
    if (SetEmptyInfo("Z", xyzabg_[Box::Z], "X", xyzabg_[Box::X], "Y", xyzabg_[Box::Y])) return 1;
  }
  if (xyzabg_[Box::ALPHA] <= 0) {
    if (SetEmptyInfo("alpha", xyzabg_[Box::ALPHA], "beta", xyzabg_[Box::BETA], "gamma", xyzabg_[Box::GAMMA])) return 1;
  }
  if (xyzabg_[Box::BETA] <= 0) {
    if (SetEmptyInfo("beta", xyzabg_[Box::BETA], "alpha", xyzabg_[Box::ALPHA], "gamma", xyzabg_[Box::GAMMA])) return 1;
  }
if (xyzabg_[Box::GAMMA] <= 0) {
    if (SetEmptyInfo("gamma", xyzabg_[Box::GAMMA], "alpha", xyzabg_[Box::ALPHA], "beta", xyzabg_[Box::BETA])) return 1;
  }

  return 0;
}

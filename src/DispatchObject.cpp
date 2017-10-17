#include "DispatchObject.h"

/** 0 are hidden categories (i.e. should not appear in help). */
const char* DispatchObject::ObjKeyword(Otype typeIn) {
  switch (typeIn) {
    case NONE: return 0;
    case PARM: return "Topology";
    case TRAJ: return "Trajectory";
    case COORDS: return "Coords";
    case ACTION: return "Action";
    case ANALYSIS: return "Analysis";
    case GENERAL: return "General";
    case SYSTEM: return "System";
    case CONTROL: return "Control";
    case HIDDEN: return 0;
    case DEPRECATED: return 0;
  }
  return 0;
}

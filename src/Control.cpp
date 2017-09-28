#include "Control.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

void Control_For_Mask::Help() const {
  mprintf("\t{atoms|residues|molecules} <var> in <mask> %s\n", DataSetList::TopIdxArgs);
}

int Control_For_Mask::SetupControl(CpptrajState& State, ArgList& argIn) {
  mask_.ResetMask();
  varname_.clear();
  commands_.clear();
  varType_ = UNKNOWN;
  // formask {atoms|residues|molecules} <var> inmask <mask> [TOP KEYWORDS]
  std::string inmask_arg = argIn.GetStringKey("in");
  if (!inmask_arg.empty()) {
    if (mask_.SetMaskString( inmask_arg )) return 1;
    if (argIn.hasKey("atoms")) varType_ = ATOMS;
    else if (argIn.hasKey("residues")) varType_ = RESIDUES;
    else if (argIn.hasKey("molecules")) varType_ = MOLECULES;
    if (varType_ == UNKNOWN) {
      mprinterr("Error: One of {atoms|residues|molecules} not specfied.\n");
      return 1;
    }
    Topology* top = State.DSL().GetTopByIndex( argIn );
    if (top == 0) return 1;
    varname_ = argIn.GetStringNext();
    if (varname_.empty()) {
      mprinterr("Error: 'for inmask': missing variable name.\n");
      return 1;
    }
    varname_ = "$" + varname_;
    if (top->SetupIntegerMask( mask_ )) return 1;
    mask_.MaskInfo();
    if (mask_.None()) return 1;
    static const char* TypeStr[] = { "ATOMS ", "RESIDUES ", "MOLECULES " };
    description_.assign("for " + std::string(TypeStr[varType_]) +
                        varname_ + " inmask " + mask_.MaskExpression() + " do");
  }
  if (varname_.empty()) {
    mprinterr("Error: for: No variable name/loop type specified.\n");
    return 1;
  }
  return 0;
}

void Control_For_Mask::Start(Varray& CurrentVars) {
  atom_ = mask_.begin();
  // Init CurrentVars
  CurrentVars.push_back( VarPair(varname_, "") );
  mprintf("DEBUG: Start: CurrentVars:");
  for (Varray::const_iterator vp = CurrentVars.begin(); vp != CurrentVars.end(); ++vp)
    mprintf(" %s=%s", vp->first.c_str(), vp->second.c_str());
  mprintf("\n");
}

Control::DoneType Control_For_Mask::CheckDone(Varray& CurrentVars) {
  if (atom_ == mask_.end()) return DONE;
  std::string atomStr = "@" + integerToString(*atom_ + 1);
  mprintf("DEBUG: Control_For_Mask: %s\n", atomStr.c_str());
  // Update CurrentVars
  Varray::iterator it = CurrentVars.begin();
  for (; it != CurrentVars.end(); ++it) {
    if (it->first == varname_) {
      it->second = atomStr;
      break;
    }
  }
  ++atom_;
  return NOT_DONE;
}

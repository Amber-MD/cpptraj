#include "Control.h"
#include "CpptrajStdio.h"

void Control_For::Help() const {
  mprintf("\t{atoms|residues|molecules} <var> inmask <mask> %s\n", DataSetList::TopIdxArgs);
}

int Control_For::SetupControl(CpptrajState& State, ArgList& argIn) {
  mask_.ResetMask();
  varname_.clear();
  commands_.clear();
  varType_ = UNKNOWN;
  // for {atoms|residues|molecules} <var> inmask <mask> [TOP KEYWORDS]
  std::string inmask_arg = argIn.GetStringKey("inmask");
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

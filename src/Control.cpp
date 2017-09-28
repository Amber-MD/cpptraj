#include "Control.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

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

void Control_For::Start() {
  atom_ = mask_.begin();
}

Control::DoneType Control_For::CheckDone() {
  if (atom_ == mask_.end()) return DONE;
  std::string atomStr = "@" + integerToString(*atom_ + 1);
  mprintf("DEBUG: Control_For: %s\n", atomStr.c_str());
  // Replace varname_ in commands with atom
  modified_commands_.clear();
  for (const_iterator it = commands_.begin(); it != commands_.end(); ++it)
  {
    modified_commands_.push_back( *it );
    for (int i = 0; i < modified_commands_.back().Nargs(); i++)
    {
      if (modified_commands_.back()[i][0] == '$') {
        if (modified_commands_.back()[i] == varname_) {
          modified_commands_.back().ChangeArg(i, atomStr);
        } else {
          mprinterr("Error: Unrecognized variable in command: %s\n",
                    modified_commands_.back()[i].c_str());
          return ERROR;
        }
      }
    }
  }
  ++atom_;
  return NOT_DONE;
}

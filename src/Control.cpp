#include "Control.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

void Control_For_Mask::Help() const {
  mprintf("\t{atoms|residues|molecules} <var> in <mask> %s\n", DataSetList::TopIdxArgs);
}

int Control_For_Mask::SetupControl(CpptrajState& State, ArgList& argIn) {
  Masks_.clear();
  description_.clear();
  Topology* currentTop = 0;
  // formask {atoms|residues|molecules} <var> inmask <mask> [TOP KEYWORDS] ...
  std::string inmask_arg = argIn.GetStringKey("in");
  while (!inmask_arg.empty()) {
    Masks_.push_back( MaskHolder() );
    MaskHolder& MH = Masks_.back();
    if (MH.mask_.SetMaskString( inmask_arg )) return 1;
    if (argIn.hasKey("atoms")) MH.varType_ = ATOMS;
    else if (argIn.hasKey("residues")) MH.varType_ = RESIDUES;
    else if (argIn.hasKey("molecules")) MH.varType_ = MOLECULES;
    if (MH.varType_ == UNKNOWN) {
      mprinterr("Error: One of {atoms|residues|molecules} not specfied.\n");
      return 1;
    }
    Topology* top = State.DSL().GetTopByIndex( argIn );
    if (top != 0) currentTop = top;
    if (currentTop == 0) return 1;
    MH.varname_ = argIn.GetStringNext();
    if (MH.varname_.empty()) {
      mprinterr("Error: 'for inmask': missing variable name.\n");
      return 1;
    }
    MH.varname_ = "$" + MH.varname_;
    if (currentTop->SetupIntegerMask( MH.mask_ )) return 1;
    MH.mask_.MaskInfo();
    if (MH.mask_.None()) return 1;
    static const char* TypeStr[] = { "ATOMS ", "RESIDUES ", "MOLECULES " };
    description_.append("for " + std::string(TypeStr[MH.varType_]) +
                        MH.varname_ + " inmask " + MH.mask_.MaskExpression() + " do ");
    inmask_arg = argIn.GetStringKey("in");
  }

  return 0;
}

void Control_For_Mask::Start(Varray& CurrentVars) {
  for (Marray::iterator MH = Masks_.begin(); MH != Masks_.end(); ++MH) {
    MH->atom_ = MH->mask_.begin();
    // Init CurrentVars
    CurrentVars.push_back( VarPair(MH->varname_, "") );
  }
  mprintf("DEBUG: Start: CurrentVars:");
  for (Varray::const_iterator vp = CurrentVars.begin(); vp != CurrentVars.end(); ++vp)
    mprintf(" %s=%s", vp->first.c_str(), vp->second.c_str());
  mprintf("\n");
}

Control::DoneType Control_For_Mask::CheckDone(Varray& CurrentVars) {
  for (Marray::iterator MH = Masks_.begin(); MH != Masks_.end(); ++MH) {
    // Exit as soon as one is done TODO check all
    if (MH->atom_ == MH->mask_.end()) return DONE;
    std::string atomStr = "@" + integerToString(*(MH->atom_) + 1);
    mprintf("DEBUG: Control_For_Mask: %s\n", atomStr.c_str());
    // Update CurrentVars
    Varray::iterator it = CurrentVars.begin();
    for (; it != CurrentVars.end(); ++it) {
      if (it->first == MH->varname_) {
        it->second = atomStr;
        break;
      }
    }
    ++(MH->atom_);
  }
  return NOT_DONE;
}

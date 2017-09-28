#include <algorithm> // std::min
#include "Control.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

void Control_For_Mask::Help() const {
  mprintf("\t{atoms|residues|molecules} <var> in <mask> %s\n", DataSetList::TopIdxArgs);
}

/** Set up each mask. */
int Control_For_Mask::SetupControl(CpptrajState& State, ArgList& argIn) {
  Masks_.clear();
  Topology* currentTop = 0;
  static const char* TypeStr[] = { "ATOMS ", "RESIDUES ", "MOLECULES " };
  description_.assign("for ");
  // formask {atoms|residues|molecules} <var> inmask <mask> [TOP KEYWORDS] ...
  std::string inMaskArg = argIn.GetStringKey("in");
  int Niterations = -1;
  while (!inMaskArg.empty()) {
    Masks_.push_back( MaskHolder() );
    MaskHolder& MH = Masks_.back();
    AtomMask currentMask;
    if (currentMask.SetMaskString( inMaskArg )) return 1;
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
    if (currentTop->SetupIntegerMask( currentMask )) return 1;
    currentMask.MaskInfo();
    if (currentMask.None()) return 1;
    // Set up indices
    if (MH.varType_ == ATOMS)
      MH.Idxs_ = currentMask.Selected();
    else if (MH.varType_ == RESIDUES) {
      int curRes = -1;
      for (AtomMask::const_iterator at = currentMask.begin(); at != currentMask.end(); ++at) {
        int res = (*currentTop)[*at].ResNum();
        if (res != curRes) {
          MH.Idxs_.push_back( res );
          curRes = res;
        }
      }
    } else if (MH.varType_ == MOLECULES) {
      int curMol = -1;
      for (AtomMask::const_iterator at = currentMask.begin(); at != currentMask.end(); ++at) {
        int mol = (*currentTop)[*at].MolNum();
        if (mol != curMol) {
          MH.Idxs_.push_back( mol );
          curMol = mol;
        }
      }
    }
    // Check number of values
    if (Niterations == -1)
      Niterations = (int)MH.Idxs_.size();
    else if ((int)MH.Idxs_.size() != Niterations)
      mprintf("Warning: # iterations %zu != previous # iterations %i\n",
              MH.Idxs_.size(), Niterations);
    if (description_ != "for ") description_.append(", ");
    description_.append(std::string(TypeStr[MH.varType_]) +
                        MH.varname_ + " in " + currentMask.MaskExpression());
    inMaskArg = argIn.GetStringKey("in");
  }
  description_.append(" do");

  return 0;
}

/** For each mask add variable to CurrentVars and initialize iterator. */
void Control_For_Mask::Start(Varray& CurrentVars) {
  for (Marray::iterator MH = Masks_.begin(); MH != Masks_.end(); ++MH) {
    MH->idx_ = MH->Idxs_.begin();
    // Init CurrentVars
    CurrentVars.push_back( VarPair(MH->varname_, "") );
  }
}

/** For each mask check if done, then update CurrentVars, then increment. */
Control::DoneType Control_For_Mask::CheckDone(Varray& CurrentVars) {
  for (Marray::iterator MH = Masks_.begin(); MH != Masks_.end(); ++MH) {
    // Exit as soon as one is done TODO check all?
    if (MH->idx_ == MH->Idxs_.end()) return DONE;
    // Get variable value
    static const char* prefix[] = {"@", ":", "^"};
    std::string maskStr = prefix[MH->varType_] + integerToString(*(MH->idx_) + 1);
    mprintf("DEBUG: Control_For_Mask: %s\n", maskStr.c_str());
    // Update CurrentVars
    Varray::iterator it = CurrentVars.begin();
    for (; it != CurrentVars.end(); ++it) {
      if (it->first == MH->varname_) {
        it->second = maskStr;
        break;
      }
    }
    ++(MH->idx_);
  }
  return NOT_DONE;
}

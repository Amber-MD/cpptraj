#include <algorithm> // std::min
#include "Control.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

void Control_For::Help() const {
  mprintf("\t{atoms|residues|molecules} <var> in <mask> %s\n", DataSetList::TopIdxArgs);
}

/** Set up each mask. */
int Control_For::SetupControl(CpptrajState& State, ArgList& argIn) {
  mprintf("    Setting up 'for' loop.\n");
  Masks_.clear();
  Topology* currentTop = 0;
  static const char* TypeStr[] = { "ATOMS ", "RESIDUES ", "MOLECULES " };
  description_.assign("for (");
  //     [<var>=<start>;<var><OP><end>;<var><OP>[<value>]]
  int Niterations = -1;
  int iarg = 0;
  while (iarg < argIn.Nargs())
  {
    // Advance to next unmarked argument.
    while (argIn.Marked(iarg) && iarg < argIn.Nargs()) iarg++;
    if (iarg == argIn.Nargs()) break;
    // Determine 'for' type
    ForType ftype = UNKNOWN;
    if      ( argIn[iarg] == "atoms"     ) ftype = ATOMS;
    else if ( argIn[iarg] == "residues"  ) ftype = RESIDUES;
    else if ( argIn[iarg] == "molecules" ) ftype = MOLECULES;
    if (ftype == UNKNOWN) {
      mprinterr("Error: One of {atoms|residues|molecules} not specfied.\n");
      return 1;
    }
    argIn.MarkArg(iarg);
    // Set up for specific type
    if (ftype == ATOMS || ftype == RESIDUES || ftype == MOLECULES)
    {
      // {atoms|residues|molecules} <var> inmask <mask> [TOP KEYWORDS]
      if (argIn[iarg+2] != "inmask") {
        mprinterr("Error: Expected 'inmask', got %s\n", argIn[iarg+2].c_str());
        return 1;
      }
      Masks_.push_back( MaskHolder() );
      MaskHolder& MH = Masks_.back();
      AtomMask currentMask;
      if (currentMask.SetMaskString( argIn.GetStringKey("inmask") )) return 1;
      MH.varType_ = ftype;
      Topology* top = State.DSL().GetTopByIndex( argIn );
      if (top != 0) currentTop = top;
      if (currentTop == 0) return 1;
      MH.varname_ = argIn.GetStringNext();
      if (MH.varname_.empty()) {
        mprinterr("Error: 'for inmask': missing variable name.\n");
        return 1;
      }
      MH.varname_ = "$" + MH.varname_;
      // Set up mask
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
      else if ((int)MH.Idxs_.size() != Niterations) {
        mprintf("Warning: # iterations %zu != previous # iterations %i\n",
                MH.Idxs_.size(), Niterations);
        Niterations = std::min((int)MH.Idxs_.size(), Niterations);
      }
      if (description_ != "for (") description_.append(", ");
      description_.append(std::string(TypeStr[MH.varType_]) +
                        MH.varname_ + " inmask " + currentMask.MaskExpression());
    } // END for mask setup
  }
  mprintf("\tLoop will execute for %i iterations.\n", Niterations);
  description_.append(") do");

  return 0;
}

/** For each mask add variable to CurrentVars and initialize iterator. */
void Control_For::Start(Varray& CurrentVars) {
  for (Marray::iterator MH = Masks_.begin(); MH != Masks_.end(); ++MH) {
    MH->idx_ = MH->Idxs_.begin();
    // Init CurrentVars
    CurrentVars.push_back( VarPair(MH->varname_, "") );
  }
}

/** For each mask check if done, then update CurrentVars, then increment. */
Control::DoneType Control_For::CheckDone(Varray& CurrentVars) {
  for (Marray::iterator MH = Masks_.begin(); MH != Masks_.end(); ++MH) {
    // Exit as soon as one is done TODO check all?
    if (MH->idx_ == MH->Idxs_.end()) return DONE;
    // Get variable value
    static const char* prefix[] = {"@", ":", "^"};
    std::string maskStr = prefix[MH->varType_] + integerToString(*(MH->idx_) + 1);
    //mprintf("DEBUG: Control_For: %s\n", maskStr.c_str());
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

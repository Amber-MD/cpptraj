#include "ForLoop_mask.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
#include "CpptrajState.h"
#include "StringRoutines.h"
#include "DataSetList.h"

/** CONSTRUCTOR */
ForLoop_mask::ForLoop_mask() :
  mtype_(NTYPES),
  currentTop_(0)
{}

void ForLoop_mask::helpText() {
  mprintf("\t{atoms|residues|molecules|molfirstres|mollastres}\n"
          "\t  <var> inmask <mask> [%s]\n", DataSetList::TopIdxArgs);
  mprintf("  Loop over elements selected by specified mask expression.\n");
}

int ForLoop_mask::SetupFor(CpptrajState& State, ArgList& argIn) {
  static const char* TypeStr[NTYPES] = { "ATOMS ", "RESIDUES ", "MOLECULES ",
                                   "MOL_FIRST_RES ", "MOL_LAST_RES " };
  // {atoms|residues|molecules} <var> inmask <mask> [TOP KEYWORDS]
  currentTop_ = 0;
  mtype_ = NTYPES;
  if      (argIn.hasKey("atoms")) mtype_ = ATOMS;
  else if (argIn.hasKey("residues")) mtype_ = RESIDUES;
  else if (argIn.hasKey("molecules")) mtype_ = MOLECULES;
  else if (argIn.hasKey("molfirstres")) mtype_ = MOLFIRSTRES;
  else if (argIn.hasKey("mollastres")) mtype_ = MOLLASTRES;
  if (mtype_ == NTYPES) {
    mprinterr("Error: No recognized type for mask 'for' loop.\n");
    return 1;
  }
  //if (argIn[iarg+2] != "inmask") {
  //  mprinterr("Error: Expected 'inmask', got %s\n", argIn[iarg+2].c_str());
  //  return 1;
  //}
  if (!argIn.Contains("inmask")) {
    mprinterr("Error: mask for loop does not contain 'inmask' keyword.\n");
    return 1;
  }
  maskExpr_ = argIn.GetStringKey("inmask");
  if (maskExpr_.empty()) {
    mprinterr("Error: Must set mask expression via 'inmask'.\n");
    return 1;
  }
  Topology* top = State.DSL().GetTopByIndex( argIn );
  if (top != 0) currentTop_ = top;
  if (currentTop_ == 0) return 1;
  // Get the variable name
  if (SetupLoopVar( State.DSL(), argIn.GetStringNext() )) return 1;

  std::string description("(" + std::string(TypeStr[mtype_]) +
                          VarName() + " inmask " + maskExpr_ + ")");
  SetDescription( description );
  return 0;
}

/** Set up the atom mask, obtain indices. */
int ForLoop_mask::BeginFor(DataSetList const& DSL) {
  if (currentTop_ == 0) {
    mprinterr("Internal Error: ForLoop_mask::BeginFor() called before Setup().\n");
    return LOOP_ERROR;
  }
  Idxs_.clear();
  AtomMask currentMask;
  // Do variable replacement
  std::string maskExpr2;
  int nReplaced = DSL.ReplaceVariables( maskExpr2, maskExpr_ );
  if (nReplaced > 0) {
    mprintf("DEBUG: old mask '%s' new mask '%s'\n", maskExpr_.c_str(), maskExpr2.c_str());
    if (currentMask.SetMaskString( maskExpr2 )) return LOOP_ERROR;
  } else {
    if (currentMask.SetMaskString( maskExpr_ )) return LOOP_ERROR;
  }
  // Set up mask
  if (currentTop_->SetupIntegerMask( currentMask )) return LOOP_ERROR;
  currentMask.MaskInfo();
  if (currentMask.None()) return 0; // No iterations
  // Set up indices
  if (mtype_ == ATOMS)
    Idxs_ = currentMask.Selected();
  else if (mtype_ == RESIDUES) {
    int curRes = -1;
    for (AtomMask::const_iterator at = currentMask.begin(); at != currentMask.end(); ++at) {
      int res = (*currentTop_)[*at].ResNum();
      if (res != curRes) {
        Idxs_.push_back( res );
        curRes = res;
      }
    }
  } else if (mtype_ == MOLECULES ||
             mtype_ == MOLFIRSTRES ||
             mtype_ == MOLLASTRES)
  {
    int curMol = -1;
    for (AtomMask::const_iterator at = currentMask.begin(); at != currentMask.end(); ++at) {
      int mol = (*currentTop_)[*at].MolNum();
      if (mol != curMol) {
        if (mtype_ == MOLECULES)
          Idxs_.push_back( mol );
        else {
          int res;
          if (mtype_ == MOLFIRSTRES)
            res = (*currentTop_)[ currentTop_->Mol( mol ).MolUnit().Front() ].ResNum();
          else // MOLLASTRES
            res = (*currentTop_)[ currentTop_->Mol( mol ).MolUnit().Back()-1 ].ResNum();
          Idxs_.push_back( res );
        }
        curMol = mol;
      }
    }
  }

  idx_ = Idxs_.begin();
  return (int)Idxs_.size();
}

bool ForLoop_mask::EndFor(DataSetList& DSL) {
  static const char* prefix[NTYPES] = {"@", ":", "^", ":", ":"};
  if (idx_ == Idxs_.end()) return true;
  // Get variable value
  std::string maskStr = prefix[mtype_] + integerToString(*(idx_) + 1);
  //mprintf("DEBUG: ControlBlock_For: %s\n", maskStr.c_str());
  // Update CurrentVars
  DSL.UpdateStringVar( VarName(), maskStr );
  // Increment
  ++(idx_);
  return false;
}

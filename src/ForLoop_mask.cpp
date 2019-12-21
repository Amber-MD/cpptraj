#include "ForLoop_mask.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
#include "CpptrajState.h"
#include "StringRoutines.h"
#include "DataSetList.h"

int ForLoop_mask::SetupFor(CpptrajState& State, std::string const& expr, ArgList& argIn) {
  static const char* TypeStr[NTYPES] = { "ATOMS ", "RESIDUES ", "MOLECULES ",
                                   "MOL_FIRST_RES ", "MOL_LAST_RES " };
  // {atoms|residues|molecules} <var> inmask <mask> [TOP KEYWORDS]
  Topology* currentTop = 0;
  mtype_ = NTYPES;
  if      ( expr == "atoms"       ) mtype_ = ATOMS;
  else if ( expr == "residues"    ) mtype_ = RESIDUES;
  else if ( expr == "molecules"   ) mtype_ = MOLECULES;
  else if ( expr == "molfirstres" ) mtype_ = MOLFIRSTRES;
  else if ( expr == "mollastres"  ) mtype_ = MOLLASTRES;
  else {
    mprinterr("Error: Unrecognized mask for loop type: %s\n", expr.c_str());
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
  AtomMask currentMask;
  if (currentMask.SetMaskString( argIn.GetStringKey("inmask") )) return 1;
  Topology* top = State.DSL().GetTopByIndex( argIn );
  if (top != 0) currentTop = top;
  if (currentTop == 0) return 1;
  SetupLoopVar( State.DSL(), argIn.GetStringNext() );
  // Set up mask
  if (currentTop->SetupIntegerMask( currentMask )) return 1;
  currentMask.MaskInfo();
  if (currentMask.None()) return 1;
  // Set up indices
  if (mtype_ == ATOMS)
    Idxs_ = currentMask.Selected();
  else if (mtype_ == RESIDUES) {
    int curRes = -1;
    for (AtomMask::const_iterator at = currentMask.begin(); at != currentMask.end(); ++at) {
      int res = (*currentTop)[*at].ResNum();
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
      int mol = (*currentTop)[*at].MolNum();
      if (mol != curMol) {
        if (mtype_ == MOLECULES)
          Idxs_.push_back( mol );
        else {
          int res;
          if (mtype_ == MOLFIRSTRES)
            res = (*currentTop)[ currentTop->Mol( mol ).BeginAtom() ].ResNum();
          else // MOLLASTRES
            res = (*currentTop)[ currentTop->Mol( mol ).EndAtom()-1 ].ResNum();
          Idxs_.push_back( res );
        }
        curMol = mol;
      }
    }
  }
  std::string description(std::string(TypeStr[mtype_]) +
                          VarName() + " inmask " + currentMask.MaskExpression());
  SetDescription( description );
  return 0;
}

int ForLoop_mask::BeginFor(DataSetList const& DSL) {
  idx_ = Idxs_.begin();
  return (int)Idxs_.size();
}

bool ForLoop_mask::EndFor(DataSetList const& DSL) {
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

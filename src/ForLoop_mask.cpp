#include "ForLoop_mask.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
#include "CpptrajState.h"

int ForLoop_mask::SetupFor(CpptrajState& State, std::string const& expr, ArgList& argIn) {
  int Niterations = -1;
  static const char* TypeStr[NTYPES] = { "ATOMS ", "RESIDUES ", "MOLECULES ",
                                   "MOL_FIRST_RES ", "MOL_LAST_RES " };
  // {atoms|residues|molecules} <var> inmask <mask> [TOP KEYWORDS]
  Topology* currentTop = 0;
  MaskType mtype = NTYPES;
  if      ( expr == "atoms"       ) mtype = ATOMS;
  else if ( expr == "residues"    ) mtype = RESIDUES;
  else if ( expr == "molecules"   ) mtype = MOLECULES;
  else if ( expr == "molfirstres" ) mtype = MOLFIRSTRES;
  else if ( expr == "mollastres"  ) mtype = MOLLASTRES;
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
  SetVarName( argIn.GetStringNext() );
  // Set up mask
  if (currentTop->SetupIntegerMask( currentMask )) return 1;
  currentMask.MaskInfo();
  if (currentMask.None()) return 1;
  // Set up indices
  if (mtype == ATOMS)
    Idxs_ = currentMask.Selected();
  else if (mtype == RESIDUES) {
    int curRes = -1;
    for (AtomMask::const_iterator at = currentMask.begin(); at != currentMask.end(); ++at) {
      int res = (*currentTop)[*at].ResNum();
      if (res != curRes) {
        Idxs_.push_back( res );
        curRes = res;
      }
    }
  } else if (mtype == MOLECULES ||
             mtype == MOLFIRSTRES ||
             mtype == MOLLASTRES)
  {
    int curMol = -1;
    for (AtomMask::const_iterator at = currentMask.begin(); at != currentMask.end(); ++at) {
      int mol = (*currentTop)[*at].MolNum();
      if (mol != curMol) {
        if (mtype == MOLECULES)
          Idxs_.push_back( mol );
        else {
          int res;
          if (mtype == MOLFIRSTRES)
            res = (*currentTop)[ currentTop->Mol( mol ).BeginAtom() ].ResNum();
          else // MOLLASTRES
            res = (*currentTop)[ currentTop->Mol( mol ).EndAtom()-1 ].ResNum();
          Idxs_.push_back( res );
        }
        curMol = mol;
      }
    }
  }
  Niterations = (int)Idxs_.size();
  std::string description(std::string(TypeStr[mtype]) +
                          VarName() + " inmask " + currentMask.MaskExpression());
  SetDescription( description );
  SetNiterations( Niterations );
  return 0;
}

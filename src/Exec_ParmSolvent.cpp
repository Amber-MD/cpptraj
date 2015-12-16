#include "Exec_ParmSolvent.h"
#include "CpptrajStdio.h"

void Exec_ParmSolvent::Help() const {
  mprintf("\t[%s] { <mask> | none }\n", DataSetList::TopIdxArgs);
  mprintf("  Set solvent for the specified topology (default first) based on <mask>.\n"
          "  If 'none' specified, remove all solvent information.\n");
}

Exec::RetType Exec_ParmSolvent::Execute(CpptrajState& State, ArgList& argIn) {
  std::string maskexpr;
  if (!argIn.hasKey("none")) {
    maskexpr = argIn.GetMaskNext();
    if ( maskexpr.empty() ) {
      mprinterr("Error: solvent: No mask specified.\n");
      return CpptrajState::ERR;
    }
  }
  // Get parm index
  Topology* parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  parm->SetSolvent( maskexpr );
  return CpptrajState::OK;
}

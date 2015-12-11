#include "Exec_ScaleDihedralK.h"
#include "CpptrajStdio.h"

void Exec_ScaleDihedralK::Help() const {
  mprintf("\t[%s] <scale factor> [<mask> [useall]]\n", DataSetList::TopArgs);
}

Exec::RetType Exec_ScaleDihedralK::Execute(CpptrajState& State, ArgList& argIn) {
  Topology* parm = State.DSL()->GetTopology( argIn );
  if (parm == 0) {
    mprinterr("Error: No topologies loaded.\n");
    return CpptrajState::ERR;
  }
  double scale_factor = argIn.getNextDouble(1.0);
  std::string maskexpr = argIn.GetMaskNext();
  bool useAll = argIn.hasKey("useall");
  mprintf("\tScaling dihedral force constants in %s by %f\n", parm->c_str(), scale_factor);
  if (!maskexpr.empty()) {
    if (useAll)
      mprintf("\tAll atoms in mask '%s' must be present to select dihedral.\n",maskexpr.c_str());
    else
      mprintf("\tAny atom in mask '%s' will select a dihedral.\n",maskexpr.c_str());
  }
  parm->ScaleDihedralK( scale_factor, maskexpr, useAll );
  return CpptrajState::OK;
}

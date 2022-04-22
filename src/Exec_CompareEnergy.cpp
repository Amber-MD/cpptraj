#include "Exec_CompareEnergy.h"
#include "CpptrajStdio.h"

// Exec_CompareEnergy::Help()
void Exec_CompareEnergy::Help() const
{

}

DataSet_Coords* Exec_CompareEnergy::GetCoordsSet(DataSetList const& DSL,
                                                 std::string const& setname)
{
  if (setname.empty()) {
    mprinterr("Error: crdout: Specify COORDS dataset name.\n");
    return 0;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)DSL.FindSetOfGroup( setname, DataSet::COORDINATES );
  if (CRD == 0) {
    mprinterr("Error: crdout: No COORDS set with name %s found.\n", setname.c_str());
    return 0;
  }
  mprintf("\tUsing set '%s'\n", CRD->legend());
  return CRD;
}

// Exec_CompareEnergy::Execute()
Exec::RetType Exec_CompareEnergy::Execute(CpptrajState& State, ArgList& argIn)
{

}

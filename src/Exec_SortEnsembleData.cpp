#include "Exec_SortEnsembleData.h"
#include "CpptrajStdio.h"

// Exec_SortEnsembleData::Help()
void Exec_SortEnsembleData::Help() const
{

}

// Exec_SortEnsembleData::Execute()
Exec::RetType Exec_SortEnsembleData::Execute(CpptrajState& State, ArgList& argIn)
{
  rprintf("DEBUG: Entering sortensembledata.\n");
  DataSetList setsToSort;
  std::string dsarg = argIn.GetStringNext();
  while (!dsarg.empty()) {
    setsToSort += State.DSL().GetMultipleSets( dsarg );
    dsarg = argIn.GetStringNext();
  }
  setsToSort.List();
  // TODO error check

  

  return CpptrajState::OK;
}

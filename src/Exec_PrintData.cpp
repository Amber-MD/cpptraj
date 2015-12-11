#include "Exec_PrintData.h"
#include "CpptrajStdio.h"

void Exec_PrintData::Help() const {
  mprintf("\t<data set>\n"
          "  Print data from data set to screen.\n");
}

Exec::RetType Exec_PrintData::Execute(CpptrajState& State, ArgList& argIn) {
  DataFile ToStdout;
  ToStdout.SetupStdout(argIn, State.Debug());
  DataSetList selected;
  std::string ds_arg = argIn.GetStringNext();
  while (!ds_arg.empty()) {
    selected += State.DSL()->GetMultipleSets( ds_arg );
    ds_arg = argIn.GetStringNext();
  }
  for (DataSetList::const_iterator ds = selected.begin(); ds != selected.end(); ++ds)
    ToStdout.AddDataSet( *ds );
  ToStdout.WriteDataOut();
  return CpptrajState::OK;
}

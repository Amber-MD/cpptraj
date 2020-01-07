#include "Exec_Show.h"
#include "CpptrajStdio.h"
#include "DataSet_StringVar.h"

void Exec_Show::Help() const {
  mprintf("  Show all current script variables and their values.\n");
}

Exec::RetType Exec_Show::Execute(CpptrajState& State, ArgList& argIn)
{
  std::string varname = argIn.GetStringNext();
  if (varname.empty()) {
    State.DSL().ListStringVar();
  } else {
    DataSetList Svars = State.DSL().GetMultipleSets( varname );
    if (Svars.empty()) {
      return CpptrajState::ERR;
    }
    for (DataSetList::const_iterator ds = Svars.begin(); ds != Svars.end(); ++ds) {
      if ( (*ds)->Type() == DataSet::STRINGVAR ) {
        DataSet_StringVar const& var = static_cast<DataSet_StringVar const&>( *(*ds) );
        mprintf("\t%s = '%s'\n", var.legend(), var.Value().c_str());
      }
    }
  }
  return CpptrajState::OK;
}

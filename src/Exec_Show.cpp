#include "Exec_Show.h"
#include "CpptrajStdio.h"
#include "DataSet_StringVar.h"

void Exec_Show::Help() const {
  mprintf("\t[<var1> ...]\n"
          "  If no variable names specified, show all current script variables and\n"
          "  their values. Otherwise, show the values of the specified script\n"
          "  variables.\n");
}

Exec::RetType Exec_Show::Execute(CpptrajState& State, ArgList& argIn)
{
  std::string varname = argIn.GetStringNext();
  if (varname.empty()) {
    State.DSL().ListStringVar();
  } else {
    while (!varname.empty()) {
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
      varname = argIn.GetStringNext();
    }
  }
  return CpptrajState::OK;
}

#include "Exec_Show.h"
#include "CpptrajStdio.h"
#include "DataSet_StringVar.h"

Exec::RetType Exec_Show::Execute(CpptrajState& State, ArgList& argIn)
{
  DataSetList vars = State.DSL().GetSetsOfType("*", DataSet::STRINGVAR);
  for (DataSetList::const_iterator it = vars.begin(); it != vars.end(); ++it)
  {
    DataSet_StringVar const& var = static_cast<DataSet_StringVar const&>( *(*it) );
    mprintf("\t%s = %s\n", var.legend(), var.Value().c_str());
  }
  return CpptrajState::OK;
}

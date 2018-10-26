#include "Exec_UpdateParameters.h"
#include "CpptrajStdio.h"
#include "DataSet_Parameters.h"
#include "DataSet_Topology.h"

// Exec_UpdateParameters::Help()
void Exec_UpdateParameters::Help() const
{
  mprintf("\t%s setname <parm set>\n", DataSetList::TopArgs);
}

// Exec_UpdateParameters::Execute()
Exec::RetType Exec_UpdateParameters::Execute(CpptrajState& State, ArgList& argIn)
{
  std::string dsname = argIn.GetStringKey("setname");
  if (dsname.empty()) {
    mprinterr("Error: Specify parameter set.\n");
    return CpptrajState::ERR;
  }
  DataSet* ds = State.DSL().GetDataSet( dsname );
  if (ds == 0) {
    mprinterr("Error: Parameter data set '%s' not found.\n", dsname.c_str());
    return CpptrajState::ERR;
  }
  if (ds->Type() != DataSet::PARAMETERS && ds->Type() != DataSet::TOPOLOGY) {
    mprinterr("Error: Set '%s' is not a parameter or topology data set.\n", ds->legend());
    return CpptrajState::ERR;
  }
  Topology* dstop = State.DSL().GetTopology( argIn );
  if (dstop == 0) {
    mprinterr("Error: No topology specified.\n");
    return CpptrajState::ERR;
  }
  Topology& top = static_cast<Topology&>( *dstop );

  mprintf("\tUpdating parameters in topology '%s' using those in set '%s'\n",
          top.c_str(), ds->legend());

  if (ds->Type() == DataSet::PARAMETERS)
    top.UpdateParams(static_cast<DataSet_Parameters const&>( *ds ));
  else if (ds->Type() == DataSet::TOPOLOGY) {
    DataSet_Topology const& topds = static_cast<DataSet_Topology const&>( *ds );
    top.UpdateParams(topds.Top().GetParameters());
  } else // Sanity check
    return CpptrajState::ERR;

  return CpptrajState::OK;
}

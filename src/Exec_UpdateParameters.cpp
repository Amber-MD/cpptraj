#include "Exec_UpdateParameters.h"
#include "CpptrajStdio.h"
#include "DataSet_Parameters.h"

// Exec_UpdateParameters::Help()
void Exec_UpdateParameters::Help() const
{
  mprintf("\t%s\n", DataSetList::TopArgs);
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
  if (ds->Type() != DataSet::PARAMETERS) {
    mprinterr("Error: Set '%s' is not a parameter data set.\n", ds->legend());
    return CpptrajState::ERR;
  }
  DataSet_Parameters const& prm = static_cast<DataSet_Parameters const&>( *ds );
  Topology* dstop = State.DSL().GetTopology( argIn );
  if (dstop == 0) {
    mprinterr("Error: No topology specified.\n");
    return CpptrajState::ERR;
  }
  Topology& top = static_cast<Topology&>( *dstop );
  // Now update parameters
  ParmHolder<int> ParmIndices;
  // Bond parameters
  for (BondArray::const_iterator b = top.Bonds().begin(); b != top.Bonds().end(); ++b)
  {
    AtomTypeHolder types(2);
    types.AddName( top[b->A1()].Type() );
    types.AddName( top[b->A2()].Type() );
    ParmIndices.AddParm( types, b->Idx(), false );
  }
  for (ParmHolder<int>::const_iterator it = ParmIndices.begin(); it != ParmIndices.end(); ++it)
    mprintf("\t%s %s %i\n", *(it->first[0]), *(it->first[0]), it->second);
  return CpptrajState::OK;
}

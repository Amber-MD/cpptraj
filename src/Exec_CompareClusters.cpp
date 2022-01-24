#include "Exec_CompareClusters.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR. */
Exec_CompareClusters::Exec_CompareClusters() :
  Exec(GENERAL)
{
  SetHidden( true );
}

// Exec_CompareClusters::Help()
void Exec_CompareClusters::Help() const
{

}

/** Get cluster number vs time data set. */
DataSet* Exec_CompareClusters::getClusterSet(std::string const& dsname, DataSetList const& DSL) {
  DataSet* ds = DSL.GetDataSet( dsname );
  if (ds == 0) {
    mprinterr("Error: '%s' does not correspond to a data set.\n", dsname.c_str());
    return 0;
  }
  if (ds->Group() != DataSet::SCALAR_1D) {
    mprinterr("Error: '%s' is not a 1D scalar data set.\n", dsname.c_str());
    return 0;
  }
  return ds;
}

// Exec_CompareClusters::Execute()
Exec::RetType Exec_CompareClusters::Execute(CpptrajState& State, ArgList& argIn)
{
  DataSet* c0 = getClusterSet( argIn.GetStringKey("set"), State.DSL() );
  if (c0 == 0) return CpptrajState::ERR;
  DataSet* c1 = getClusterSet( argIn.GetStringKey("set"), State.DSL() );
  if (c1 == 0) return CpptrajState::ERR;
  mprintf("\tCluster set 0: '%s'\n", c0->legend());
  mprintf("\tCluster set 1: '%s'\n", c1->legend());

  return CpptrajState::OK;
}

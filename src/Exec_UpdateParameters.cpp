#include "Exec_UpdateParameters.h"
#include "CpptrajStdio.h"
#include "DataSet_Parameters.h"
#include "DataSet_Topology.h"
#include "Parm/AssignParams.h"
#include "Parm/GetParams.h"

const char* Exec_UpdateParameters::disclaimer_ = "Warning: This command is provided for convenience only.\nWarning: For editing topology files, ParmEd is a much better alternative.\n";

// Exec_UpdateParameters::Help()
void Exec_UpdateParameters::Help() const
{
  mprintf("\tsetname <parmset> [verbose <#>]\n"
          "\t%s\n", DataSetList::TopArgs);
  mprintf("  Update parameters in specified topology with those from <parm set>.\n"
          "  <parm set> can either be a parameter set or a topology. If a\n"
          "  parameter from <parm set> does not exist in the topology it\n"
          "  will be added.\n");
  mprintf("%s", disclaimer_);
}

// Exec_UpdateParameters::Execute()
Exec::RetType Exec_UpdateParameters::Execute(CpptrajState& State, ArgList& argIn)
{
  mprintf("%s", disclaimer_);
  int verbose = argIn.getKeyInt("verbose", 1);
  std::string dsname = argIn.GetStringKey("setname");
  //bool genAngles = argIn.hasKey("genangles");
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
  //if (genAngles)
  //  mprintf("\tWill attempt to generate angle/dihedral information from bonds.\n");

  //if (genAngles) {
  //  if (Cpptraj::Structure::GenerateBondAngleTorsionArrays( top )) {
  //    mprinterr("Error: Could not generate angle/dihedral information.\n");
  //    return CpptrajState::ERR;
  //  }
  //}

  Cpptraj::Parm::AssignParams AP;
  AP.SetDebug( State.Debug() );
  AP.SetVerbose( verbose );
  if (ds->Type() == DataSet::PARAMETERS)
    AP.UpdateParameters( top, static_cast<DataSet_Parameters const&>( *ds ) );
  else if (ds->Type() == DataSet::TOPOLOGY) {
    Cpptraj::Parm::GetParams GP;
    GP.SetDebug( State.Debug() );
    DataSet_Topology const& topds = static_cast<DataSet_Topology const&>( *ds );
    AP.UpdateParameters( top, GP.GetParameters( topds.Top() ) );
  } else // Sanity check
    return CpptrajState::ERR;

  return CpptrajState::OK;
}

#include "Exec_LoadCrd.h"
#include "CpptrajStdio.h"
#include "Trajin_Single.h"

void Exec_LoadCrd::Help() const {
  mprintf("\t<filename> [%s] [<trajin args>] [name <name>]\n", DataSetList::TopArgs);
  mprintf("  Load trajectory <filename> as a COORDS data set named <name> (default <filename>).\n");
}

Exec::RetType Exec_LoadCrd::Execute(CpptrajState& State, ArgList& argIn) {
  // Get parm
  Topology* parm = State.DSL().GetTopology( argIn );
  if (parm == 0) {
    mprinterr("Error: loadcrd: No parm files loaded.\n");
    return CpptrajState::ERR;
  }
  // Load trajectory
  Trajin_Single trajin;
  trajin.SetDebug( State.Debug() );
  if (trajin.SetupTrajRead(argIn.GetStringNext(), argIn, parm)) {
    mprinterr("Error: loadcrd: Could not set up input trajectory.\n");
    return CpptrajState::ERR;
  }
  // Create input frame
  Frame frameIn;
  frameIn.SetupFrameV(parm->Atoms(), trajin.TrajCoordInfo());
  // Set up metadata with file name and output set name
  std::string setname = argIn.GetStringKey("name");
  // NOTE: For backwards compat.
  if (setname.empty()) setname = argIn.GetStringNext();
  MetaData md( trajin.Traj().Filename(), setname, -1 );
  // Check if set already present
  DataSet_Coords* coords = 0;
  DataSet* ds = State.DSL().FindSetOfType( setname, DataSet::COORDS );
  if (ds == 0) {
    // Create Set 
    coords = (DataSet_Coords*)State.DSL().AddSet(DataSet::COORDS, md);
    if (coords == 0) {
      mprinterr("Error: loadcrd: Could not set up COORDS data set.\n");
      return CpptrajState::ERR;
    }
    coords->CoordsSetup( *parm, trajin.TrajCoordInfo() );
    mprintf("\tLoading trajectory '%s' as '%s'\n", trajin.Traj().Filename().full(),
            coords->legend());
  } else {
    // Check that set is actually coords.
    if (ds->Type() != DataSet::COORDS) {
      mprinterr("Error: Set %s present but is not of type COORDS.\n", ds->legend());
      return CpptrajState::ERR;
    }
    coords = (DataSet_Coords*)ds;
    // Check that topology matches. For now just check # atoms.
    if (parm->Natom() != coords->Top().Natom()) {
      mprinterr("Error: Trajectory '%s' # atoms %i does not match COORDS data set '%s' (%i)\n",
                trajin.Traj().Filename().full(), parm->Natom(),
                coords->legend(), coords->Top().Natom());
      return CpptrajState::ERR;
    }
    mprintf("\tAppending trajectory '%s' to COORDS data set '%s'\n",
            trajin.Traj().Filename().full(), coords->legend());
  }
  // Read trajectory TODO progress bar
  trajin.BeginTraj();
  trajin.Traj().PrintInfoLine();
  while (trajin.GetNextFrame( frameIn ))
    coords->AddFrame( frameIn );
  trajin.EndTraj();
  return CpptrajState::OK;
}

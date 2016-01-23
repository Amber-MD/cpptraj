#include "Exec_LoadTraj.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords_TRJ.h"

void Exec_LoadTraj::Help() const {
  mprintf("\tname <setname> [<filename>]\n"
          "  Create/add to TRAJ data set named <setname>. If no <filename> given, convert\n"
          "  currently loaded input trajectories to TRAJ data set; otherwise add <filename>\n"
          "  to TRAJ data set <setname>\n");
}

Exec::RetType Exec_LoadTraj::Execute(CpptrajState& State, ArgList& argIn) {
  // Get Keywords
  std::string setname = argIn.GetStringKey("name");
  if (setname.empty()) {
    mprinterr("Error: Must provide data set name ('name <setname>')\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords_TRJ* trj = (DataSet_Coords_TRJ*)
                            State.DSL().FindSetOfType(setname, DataSet::TRAJ);
  if (trj == 0)
    trj = (DataSet_Coords_TRJ*)
          State.DSL().AddSet(DataSet::TRAJ, setname, "__DTRJ__");
  if (trj == 0) {
    mprinterr("Error: Could not set up TRAJ data set.\n");
    return CpptrajState::ERR;
  }
  std::string trajname = argIn.GetStringNext();
  if (trajname.empty()) {
    // Add all existing input trajectories
    if (State.InputTrajList().empty()) {
      mprinterr("Error: No input trajectories loaded.\n");
      return CpptrajState::ERR;
    }
    if (State.Mode() != CpptrajState::NORMAL) {
      mprinterr("Error: Cannot convert ensemble input trajectories to data.\n");
      return CpptrajState::ERR;
    }
    mprintf("\tSaving currently loaded input trajectories as data set with name '%s'\n",
            setname.c_str());
    for (TrajinList::trajin_it Trajin = State.InputTrajList().trajin_begin();
                               Trajin != State.InputTrajList().trajin_end(); ++Trajin)
      if (trj->AddInputTraj( *Trajin )) return CpptrajState::ERR;
    // TODO: Clear input trajectories from trajinList?
  } else {
    // Add the named trajectory
    Topology* top = State.DSL().GetTopology(argIn);
    if (top == 0) {
      mprinterr("Error: No topologies loaded.\n");
      return CpptrajState::ERR;
    }
    if (trj->AddSingleTrajin( trajname, argIn, top ))
      return CpptrajState::ERR;
  }
  return CpptrajState::OK;
}

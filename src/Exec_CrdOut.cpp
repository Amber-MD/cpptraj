#include "Exec_CrdOut.h"
#include "Trajout_Single.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"

void Exec_CrdOut::Help() const {
  mprintf("\t<crd set> <filename> [<trajout args>] [crdframes <start>,<stop>,<offset>]\n"
          "  Write COORDS data set <crd set> to trajectory file <filename>\n");
}

Exec::RetType Exec_CrdOut::Execute(CpptrajState& State, ArgList& argIn) {
# ifdef MPI
  Exec::RetType ret = CpptrajState::OK;
  int err = 0;
  // Create a communicator that just contains the master.
  int ID = MPI_UNDEFINED;
  if (Parallel::TrajComm().Master())
    ID = Parallel::World().Rank();
  //rprintf("DEBUG: About to create new comm, ID= %i\n", ID);
  trajComm_ = Parallel::World().Split( ID );
  if (ID != MPI_UNDEFINED) {
    mprintf("Warning: '%s' command does not yet use multiple MPI processes.\n", argIn.Command());
    ret = WriteCrd(State, argIn);
    if (ret != CpptrajState::OK)
      err = 1;
  }
  trajComm_.Reset();
  if (Parallel::World().CheckError( err ))
    ret = CpptrajState::ERR;
  return ret;
# else
  return (WriteCrd(State, argIn));
# endif
}

// Exec_CrdOut::WriteCrd()
Exec::RetType Exec_CrdOut::WriteCrd(CpptrajState& State, ArgList& argIn) {
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: crdout: Specify COORDS dataset name.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL().FindSetOfGroup( setname, DataSet::COORDINATES );
  if (CRD == 0) {
    mprinterr("Error: crdout: No COORDS set with name %s found.\n", setname.c_str());
    return CpptrajState::ERR;
  }
  mprintf("\tUsing set '%s'\n", CRD->legend());
  setname = argIn.GetStringNext(); // Output traj file name
  // Start, stop, offset
  TrajFrameCounter frameCount;
  ArgList crdarg( argIn.GetStringKey("crdframes"), "," );
  if (frameCount.CheckFrameArgs( CRD->Size(), crdarg )) return CpptrajState::ERR;
  frameCount.PrintInfoLine( CRD->legend() );
  Trajout_Single outtraj;
  if (outtraj.PrepareTrajWrite( setname, argIn, State.DSL(), CRD->TopPtr(), CRD->CoordsInfo(),
                                CRD->Size(), TrajectoryFile::UNKNOWN_TRAJ))
  {
    mprinterr("Error: crdout: Could not set up output trajectory.\n");
    return CpptrajState::ERR;
  }
  outtraj.PrintInfo(0);
  Frame currentFrame = CRD->AllocateFrame();
  ProgressBar progress( frameCount.TotalReadFrames() );
  int set = 0;
  for (int frame = frameCount.Start(); frame < frameCount.Stop();
           frame += frameCount.Offset(), ++set)
  {
    progress.Update( set );
    CRD->GetFrame( frame, currentFrame );
    if ( outtraj.WriteSingle( frame, currentFrame ) ) {
      mprinterr("Error writing %s to output trajectory, frame %i.\n",
                CRD->legend(), frame + 1);
      break;
    }
  }
  return CpptrajState::OK;
}

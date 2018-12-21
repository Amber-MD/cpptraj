#include "Exec_CrdAction.h"
#include "CpptrajStdio.h"
#include "Command.h"
#include "Timer.h"
#include "ProgressBar.h"
#include "DataSet_Coords_CRD.h"

void Exec_CrdAction::Help() const {
  mprintf("\t<crd set> <actioncmd> [<action args>] [crdframes <start>,<stop>,<offset>]\n"
          "  Perform action <actioncmd> on COORDS data set <crd set>.\n");
}

Exec::RetType Exec_CrdAction::DoCrdAction(CpptrajState& State, ArgList& actionargs,
                                          DataSet_Coords* CRD, Action* act,
                                          TrajFrameCounter const& frameCount) const
{
  Timer total_time;
  total_time.Start();
# ifdef MPI
  ActionInit state(State.DSL(), State.DFL(), trajComm_);
# else
  ActionInit state(State.DSL(), State.DFL());
# endif
  if ( act->Init( actionargs, state, State.Debug() ) != Action::OK )
    return CpptrajState::ERR;
  actionargs.CheckForMoreArgs();
  // Set up frame and parm for COORDS.
  ActionSetup originalSetup( CRD->TopPtr(), CRD->CoordsInfo(), CRD->Size() );
  Frame originalFrame = CRD->AllocateFrame();
  // Set up for this topology 
  Action::RetType setup_ret = act->Setup( originalSetup );
  if ( setup_ret == Action::ERR || setup_ret == Action::SKIP )
    return CpptrajState::ERR;
  // If the topology was modified, we will need a new COORDS set.
  DataSet_Coords* crdOut = 0;
  if ( setup_ret == Action::MODIFY_TOPOLOGY ) {
    // This will not work for a TRJ set.
    switch ( CRD->Type() ) {
      case DataSet::TRAJ      : mprinterr("Error: Cannot modify TRAJ data sets.\n"); break;
      case DataSet::COORDS    : crdOut = (DataSet_Coords*)new DataSet_Coords_CRD(); break;
      case DataSet::REF_FRAME : crdOut = (DataSet_Coords*)new DataSet_Coords_REF(); break;
      default: crdOut = 0; // SANITY
    }
    if (crdOut == 0) return CpptrajState::ERR;
    mprintf("Info: crdaction: COORDS set '%s' will be modified by action '%s'\n",
            CRD->legend(), actionargs.Command());
    if (frameCount.TotalReadFrames() != (int)CRD->Size())
      mprintf("Info: crdaction: Previous size= %zu, new size is %i\n",
              CRD->Size(), frameCount.TotalReadFrames());
    // Set up set, copy original metadata
    crdOut->SetMeta( CRD->Meta() );
    if (crdOut->CoordsSetup( originalSetup.Top(), originalSetup.CoordInfo() ))
      return CpptrajState::ERR;
    DataSet::SizeArray mfArray(1, frameCount.TotalReadFrames());
    if (crdOut->Allocate( mfArray )) return CpptrajState::ERR;
  }
    
  // Loop over all frames in COORDS.
  ProgressBar* progress = 0;
  if (State.ShowProgress())
    progress = new ProgressBar( frameCount.TotalReadFrames() );
  int set = 0;
  for (int frame = frameCount.Start(); frame < frameCount.Stop();
           frame += frameCount.Offset(), ++set)
  {
    // Since Frame can be modified by actions, save original and use currentFrame
    ActionFrame frm( &originalFrame, set );
    if (progress != 0) progress->Update( set );
    CRD->GetFrame( frame, originalFrame );
    Action::RetType ret = act->DoAction( set, frm );
    if (ret == Action::ERR) {
      mprinterr("Error: crdaction: Frame %i, set %i\n", frame + 1, set + 1);
      break;
    }
    // Check if frame was modified. If so, update COORDS.
    if ( ret == Action::MODIFY_COORDS ) {
      if (crdOut != 0)
        crdOut->AddFrame( frm.Frm() );
      else
        CRD->SetCRD( frame, frm.Frm() );
    }
  } 
  if (progress != 0) delete progress;
# ifdef MPI
  act->SyncAction();
# endif
  // If topology was modified, replace old set with new.
  if ( setup_ret == Action::MODIFY_TOPOLOGY ) {
    mprintf("Info: crdaction: Topology for '%s' was modified by action '%s'\n",
            CRD->legend(), actionargs.Command());
    State.DSL().RemoveSet( CRD );
    State.DSL().AddSet( crdOut );
  } 
  act->Print();
  State.MasterDataFileWrite();
  total_time.Stop();
  mprintf("TIME: Total action execution time: %.4f seconds.\n", total_time.Total());
  return CpptrajState::OK;
}

// Exec_CrdAction::Execute()
Exec::RetType Exec_CrdAction::Execute(CpptrajState& State, ArgList& argIn) {
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
    ret = ProcessArgs(State, argIn);
    if (ret != CpptrajState::OK)
      err = 1;
  }
  trajComm_.Reset();
  if (Parallel::World().CheckError( err ))
    ret = CpptrajState::ERR;
  return ret;
# else
  return (ProcessArgs(State, argIn));
# endif
}

// Exec_CrdAction::ProcessArgs()
Exec::RetType Exec_CrdAction::ProcessArgs(CpptrajState& State, ArgList& argIn) {
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: %s: Specify COORDS dataset name.\n", argIn.Command());
    return CpptrajState::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL().FindCoordsSet( setname );
  if (CRD == 0) {
    mprinterr("Error: %s: No COORDS set with name %s found.\n", argIn.Command(), setname.c_str());
    return CpptrajState::ERR;
  }
  mprintf("\tUsing set '%s'\n", CRD->legend());
  // Start, stop, offset
  TrajFrameCounter frameCount;
  ArgList crdarg( argIn.GetStringKey("crdframes"), "," );
  if (frameCount.CheckFrameArgs( CRD->Size(), crdarg )) return CpptrajState::ERR;
  frameCount.PrintInfoLine(CRD->legend());
  ArgList actionargs = argIn.RemainingArgs();
  actionargs.MarkArg(0);
  Cmd const& cmd = Command::SearchTokenType( DispatchObject::ACTION, actionargs.Command() );
  if ( cmd.Empty() ) return CpptrajState::ERR;
  Action* act = (Action*)cmd.Alloc();
  if (act == 0) return CpptrajState::ERR;
  CpptrajState::RetType err = DoCrdAction(State, actionargs, CRD, act, frameCount);
  delete act;
  return err;
}

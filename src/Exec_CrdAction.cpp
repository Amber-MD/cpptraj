#include "Exec_CrdAction.h"
#include "CpptrajStdio.h"
#include "Command.h"
#include "Timer.h"
#include "ProgressBar.h"

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
  ActionInit state(State.DSL(), State.DFL());
  if ( act->Init( actionargs, state, State.Debug() ) != Action::OK )
    return CpptrajState::ERR;
  actionargs.CheckForMoreArgs();
  // Set up frame and parm for COORDS.
  ActionSetup originalSetup( CRD->TopPtr(), CRD->CoordsInfo(), CRD->Size() );
  Frame originalFrame = CRD->AllocateFrame();
  ActionFrame frm( &originalFrame, 0 );
  // Set up for this topology 
  Action::RetType setup_ret = act->Setup( originalSetup );
  if ( setup_ret == Action::ERR || setup_ret == Action::SKIP )
    return CpptrajState::ERR;
  // Loop over all frames in COORDS.
  ProgressBar progress( frameCount.TotalReadFrames() );
  int set = 0;
  for (int frame = frameCount.Start(); frame < frameCount.Stop();
           frame += frameCount.Offset(), ++set)
  { 
    progress.Update( set );
    CRD->GetFrame( frame, originalFrame );
    frm.SetTrajoutNum( set );
    Action::RetType ret = act->DoAction( set, frm );
    if (ret == Action::ERR) {
      mprinterr("Error: crdaction: Frame %i, set %i\n", frame + 1, set + 1);
      break;
    }
    // Check if frame was modified. If so, update COORDS.
    if ( ret == Action::MODIFY_COORDS )
      CRD->SetCRD( frame, frm.Frm() );
  }
  // Check if parm was modified. If so, update COORDS.
  if ( setup_ret == Action::MODIFY_TOPOLOGY ) {
    mprintf("Info: crdaction: Parm for %s was modified by action %s\n",
            CRD->legend(), actionargs.Command());
    CRD->CoordsSetup( originalSetup.Top(), originalSetup.CoordInfo() );
  } 
  act->Print();
  State.MasterDataFileWrite();
  total_time.Stop();
  mprintf("TIME: Total action execution time: %.4f seconds.\n", total_time.Total());
  return CpptrajState::OK;
}

Exec::RetType Exec_CrdAction::Execute(CpptrajState& State, ArgList& argIn) {
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

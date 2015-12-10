#include "Exec_CrdAction.h"
#include "CpptrajStdio.h"
#include "Command.h"

void Exec_CrdAction::Help() const {
  mprintf("\t<crd set> <actioncmd> [<action args>] [crdframes <start>,<stop>,<offset>]\n"
          "  Perform action <actioncmd> on COORDS data set <crd set>.\n");
}

Exec::RetType Exec_CrdAction::Execute(CpptrajState& State, ArgList& argIn) {
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: %s: Specify COORDS dataset name.\n", argIn.Command());
    return CpptrajState::ERR;
  }
  //DataSet_Coords* CRD = (DataSet_Coords*)State.DSL()->FindCoordsSet( setname );
  //if (CRD == 0) {
  //  mprinterr("Error: %s: No COORDS set with name %s found.\n", argIn.Command(), setname.c_str());
  //  return CpptrajState::ERR;
  //}
  //mprintf("\tUsing set '%s'\n", CRD->legend());
  // Start, stop, offset
  //TrajFrameCounter frameCount;
  ArgList crdarg( argIn.GetStringKey("crdframes"), "," );
  //if (frameCount.CheckFrameArgs( CRD->Size(), crdarg )) return CpptrajState::ERR;
  //frameCount.PrintInfoLine(CRD->legend());
  ArgList actionargs = argIn.RemainingArgs();
  actionargs.MarkArg(0);
  Cmd const& cmd = Command::SearchTokenType( DispatchObject::ACTION, actionargs.Command() );
  if ( cmd.Empty() ) return CpptrajState::ERR;
  Action* act = (Action*)cmd.Alloc();
  if (act == 0) return CpptrajState::ERR;
  delete act;
  return CpptrajState::OK;
}

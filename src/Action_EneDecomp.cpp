#include "Action_EneDecomp.h"
#include "CpptrajStdio.h"

// Action_EneDecomp::Help()
void Action_EneDecomp::Help() const {

}

// Action_EneDecomp::Init()
Action::RetType Action_EneDecomp::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  if (eneDecomp_.InitDecomposer( actionArgs, init.DSL(), debugIn ))
    return  Action::ERR;
  mprintf("    ENEDECOMP: Decomposing energy for selected atoms.\n");
  eneDecomp_.PrintOpts();

  return Action::OK;
}

// Action_EneDecomp::Setup()
Action::RetType Action_EneDecomp::Setup(ActionSetup& setup)
{
  return Action::ERR; // FIXME
}

// Action_EneDecomp::DoAction()
Action::RetType Action_EneDecomp::DoAction(int frameNum, ActionFrame& frm)
{
  return Action::ERR; // FIXME
}

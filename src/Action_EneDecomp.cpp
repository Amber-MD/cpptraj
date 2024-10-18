#include "Action_EneDecomp.h"
#include "CpptrajStdio.h"

// Action_EneDecomp::Help()
void Action_EneDecomp::Help() const {
  Cpptraj::Energy::EnergyDecomposer::HelpText(); 
}

// Action_EneDecomp::Init()
Action::RetType Action_EneDecomp::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  if (eneDecomp_.InitDecomposer( actionArgs, init.DSL(), init.DFL(), debugIn ))
    return  Action::ERR;
  mprintf("    ENEDECOMP: Decomposing energy for selected atoms.\n");
  eneDecomp_.PrintOpts();

  return Action::OK;
}

// Action_EneDecomp::Setup()
Action::RetType Action_EneDecomp::Setup(ActionSetup& setup)
{
  int ret = eneDecomp_.SetupDecomposer( setup.Top(), setup.CoordInfo().TrajBox() );
  if (ret == -1)
    return Action::SKIP;
  else if (ret == 1)
    return Action::ERR;
  return Action::OK;
}

// Action_EneDecomp::DoAction()
Action::RetType Action_EneDecomp::DoAction(int frameNum, ActionFrame& frm)
{
  if (eneDecomp_.CalcEne( frm.Frm() ))
    return Action::ERR;
  return Action::OK;
}

#ifdef MPI
int Action_EneDecomp::SyncAction() {
  return eneDecomp_.ReduceToMaster(trajComm_);
}
#endif

// Action_EneDecomp::Print()
void Action_EneDecomp::Print() {
  eneDecomp_.FinishCalc();
}

#include "Action_Keep.h"
#include "CpptrajStdio.h"
#include "DataSet_string.h"

/** CONSTRUCTOR */
Action_Keep::Action_Keep() :
  bridgeData_(0),
  nbridge_(0)
{}

// Action_Keep::Help()
void Action_Keep::Help() const {
  mprintf("\t[bridgedata <bridge data set> [nbridge <#>]]\n"
          "  Keep only specified parts of the system.\n"
         );
}

// Action_Keep::Init()
Action::RetType Action_Keep::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  std::string bridgeDataName = actionArgs.GetStringKey("bridgedata");
  if (!bridgeDataName.empty()) {
    DataSet* ds = init.DSL().GetDataSet( bridgeDataName );
    if (ds == 0) {
      mprinterr("Error: No data set corresponding to '%s'\n", bridgeDataName.c_str());
      return Action::ERR;
    }
    if (ds->Type() != DataSet::STRING) {
      mprinterr("Error: Bridge data set '%s' is not a string data set.\n", ds->legend());
      return Action::ERR;
    }
    bridgeData_ = static_cast<DataSet_string*>( ds );
    nbridge_ = actionArgs.getKeyInt("nbridge", 1);
    if (nbridge_ < 1) {
      mprinterr("Error: Number of bridging residues to keep must be >= 1.\n");
      return Action::ERR;
    }
  } else {
    mprinterr("Error: Nothing specified to keep.\n");
    return Action::ERR;
  }

  mprintf("    KEEP:\n");
  if (bridgeData_ != 0) {
    mprintf("\tBridge ID data set: %s\n", bridgeData_->legend());
    mprintf("\t# of bridging residues to keep: %i\n", nbridge_);
  }
  return Action::OK;
}

// Action_Keep::Setup()
Action::RetType Action_Keep::Setup(ActionSetup& setup)
{
  return Action::OK;
}

// Action_Keep::DoAction()
Action::RetType Action_Keep::DoAction(int frameNum, ActionFrame& frm)
{
  Action::RetType err = Action::OK;
  if (bridgeData_ != 0)
    err = keepBridge(frameNum, frm);

  return err;
}

/** Want to keep only residues specified in a bridge ID data set. */
Action::RetType Action_Keep::keepBridge(int frameNum, ActionFrame& frm) {
  // Ensure we can get data
  if ((unsigned int)frameNum >= bridgeData_->Size()) {
    mprinterr("Error: Frame # %i is out of range for bridge data '%s' (size is %zu)\n",
              frameNum+1, bridgeData_->legend(), bridgeData_->Size());
    return Action::ERR;
  }
  std::string const& bridgeIDstr = (*bridgeData_)[frameNum];
  if (bridgeIDstr == "None") {
    mprintf("DEBUG: Frame %i has no bridging waters.\n");
    return Action::SKIP;
  }
  ArgList bridgeID( bridgeIDstr, "," );
  mprintf("DEBUG: Frame %i has %i bridging waters.\n", frameNum+1, bridgeID.Nargs());

  return Action::OK;
}

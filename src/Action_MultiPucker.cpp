#include "Action_MultiPucker.h"
#include "CpptrajStdio.h"

Action_MultiPucker::Action_MultiPucker() :
  outfile_(0),
  masterDSL_(0)
{}

// Action_MultiPucker::Help()
void Action_MultiPucker::Help() const {

}

// Action_MultiPucker::Init()
Action::RetType Action_MultiPucker::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get Keywords
  outfile_ = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs);
  std::string resrange_arg = actionArgs.GetStringKey("resrange");
  if (!resrange_arg.empty())
    if (resRange_.SetRange( resrange_arg )) return Action::ERR;
  // Search for known pucker keywords
  if (puckerSearch_.SearchForArgs( actionArgs )) return Action::ERR;
  // Get custom pucker args
  if (puckerSearch_.SearchForNewTypeArgs( actionArgs )) return Action::ERR;
  // If no pucker types are yet selected, this will select all.
  puckerSearch_.SearchForAll();

  // Setup DataSet(s) name
  dsetname_ = actionArgs.GetStringNext();

  mprintf("    MULTIPUCKER: Calculating");
  puckerSearch_.PrintTypes();
  if (!resRange_.Empty())
    mprintf(" puckers for residues in range %s\n", resRange_.RangeArg());
  else
    mprintf(" puckers for all solute residues.\n");
  if (!dsetname_.empty())
    mprintf("\tDataSet name: %s\n", dsetname_.c_str());
  if (outfile_ != 0) mprintf("\tOutput to %s\n", outfile_->DataFilename().base());
  //if (minTorsion_ > -180.0) 
  //  mprintf("\tOutput range is 0 to 360 degrees.\n");
  //else
  //  mprintf("\tOutput range is -180 to 180 degrees.\n");
  init.DSL().SetDataSetsPending(true);
  masterDSL_ = init.DslPtr();
  return Action::OK;
}

// Action_MultiPucker::Setup()
Action::RetType Action_MultiPucker::Setup(ActionSetup& setup)
{

}

// Action_MultiPucker::DoAction()
Action::RetType Action_MultiPucker::DoAction(int frameNum, ActionFrame& frm)
{

}

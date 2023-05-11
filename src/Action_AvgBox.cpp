#include "Action_AvgBox.h"
#include "CpptrajStdio.h"

// Action_AvgBox::Help()
void Action_AvgBox::Help() const {
  mprintf("\t[name <setname>] [out <file>]\n"
          "  Calculate average box.\n");
}

// Action_AvgBox::Init()
Action::RetType Action_AvgBox::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  std::string dsname = actionArgs.GetStringNext();

  // Add DataSet
  boxMatrix_ = init.DSL().AddSet( DataSet::MAT3X3, MetaData(dsname, "avg") );
  if (boxMatrix_ == 0) {
    mprinterr("Error: Could not allocate data set for average box.\n");
    return Action::ERR;
  }
# ifdef MPI
  // This set does not need to be synced since averaging is done by avgbox_
  boxMatrix_->SetNeedsSync( false );
# endif
  if (outfile != 0) outfile->AddDataSet( boxMatrix_ );

  mprintf("    AVGBOX: Calculating average box.\n");
  mprintf("\tAverage box set: %s\n", boxMatrix_->legend());
  if (outfile != 0) mprintf("\tOutput to file: %s\n", outfile->DataFilename().full());
  return Action::OK;
}

// Action_AvgBox::Setup()
Action::RetType Action_AvgBox::Setup(ActionSetup& setup)
{

}

// Action_AvgBox::DoAction()
Action::RetType Action_AvgBox::DoAction(int frameNum, ActionFrame& frm)
{

}

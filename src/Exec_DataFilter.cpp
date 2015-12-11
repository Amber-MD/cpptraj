#include "Exec_DataFilter.h"
#include "CpptrajStdio.h"
#include "Action_FilterByData.h"
#include "ProgressBar.h"

void Exec_DataFilter::Help() const {
  mprintf("\t<dataset arg> min <min> max <max> [out <file> [name <setname>]]\n"
          "  Create a data set (optionally named <setname>) containing 1 for\n"
          "  data within given <min> and <max> criteria for each specified\n"
          "  data set. There must be at least one <min> and <max> argument,\n"
          "  and can be as many as there are specified data sets.\n");
}

Exec::RetType Exec_DataFilter::Execute(CpptrajState& State, ArgList& argIn) {
  Action_FilterByData filterAction;
  ActionInit state(*State.DSL(), *State.DFL());
  if (filterAction.Init(argIn, state, State.Debug()) != Action::OK)
    return CpptrajState::ERR;
  size_t nframes = filterAction.DetermineFrames();
  if (nframes < 1) {
    mprinterr("Error: No data to filter. All sets must contain some data.\n");
    return CpptrajState::ERR;
  }
  ProgressBar progress( nframes );
  ActionFrame frm;
  for (size_t frame = 0; frame != nframes; frame++) {
    progress.Update( frame );
    filterAction.DoAction(frame, frm); // Filter does not need frame.
  }
  // Trigger master datafile write just in case
  State.MasterDataFileWrite();
  return CpptrajState::OK;
}

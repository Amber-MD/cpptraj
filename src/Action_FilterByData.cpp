#include "Action_FilterByData.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
Action_FilterByData::Action_FilterByData() {}

void Action_FilterByData::Help() const {
  mprintf("\t{<dataset arg> min <min> max <max> ...} [out <file>] [name <setname>]\n"
          "\t[multi]\n"
          "  This action has two modes. In the first mode, for all following actions\n"
          "  only frames that are between <min> and <max> of all data sets selected\n"
          "  by each <dataset arg> are allowed to pass. A data set with name <setname>\n"
          "  will be created containing a 1 if the frame passed and 0 if the frame was\n"
          "  filtered out.\n"
          "  If 'multi' is specified then only filter data sets will be created for each\n"
          "  data set instead.\n"
          "  There must be at least one <min> and <max> argument, and can be as many as\n"
          "  there are specified data sets.\n");
}

// Action_FilterByData::Init()
Action::RetType Action_FilterByData::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  if (dataFilter_.InitFilter(actionArgs, init.DSL(), init.DFL(), debugIn))
    return Action::ERR;

  mprintf("    FILTER:");
  if (!dataFilter_.IsMulti())
    mprintf(" Filtering out frames using %u data sets.\n", dataFilter_.NinputSets());
  else
    mprintf(" Creating filter data sets for %zu data sets.\n", dataFilter_.NinputSets());
  dataFilter_.PrintInputSets();
  if (dataFilter_.OutputFile() != 0)
    mprintf("\tFilter frame info will be written to %s\n",
            dataFilter_.OutputFile()->DataFilename().full());
# ifdef MPI
  if (!dataFilter_.IsMulti() && init.TrajComm().Size() > 1)
    mprintf("Warning: Trajectories written after 'filter' may have issues if\n"
            "Warning:   the number of processes writing is > 1 (currently %i processes)\n",
            init.TrajComm().Size());
# endif
  return Action::OK;
}

// Action_FilterByData::DoAction()
Action::RetType Action_FilterByData::DoAction(int frameNum, ActionFrame& frm)
{
  DataFilter::ResultType result = dataFilter_.FilterIndex( frm.TrajoutNum() );

  if (result == DataFilter::FILTERED)
    return Action::SUPPRESS_COORD_OUTPUT;

  return Action::OK;
}

// Action_FilterByData::Print()
void Action_FilterByData::Print() {
  dataFilter_.Finalize();
}

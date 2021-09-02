#include "Exec_DataFilter.h"
#include "DataFilter.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"

void Exec_DataFilter::Help() const {
  DataFilter::PrintKeywords();
  mprintf("  Create a data set (optionally named <setname>) containing 1 for\n"
          "  data within given <min> and <max> criteria for each specified\n"
          "  data set. There must be at least one <min> and <max> argument,\n"
          "  and can be as many as there are specified data sets.\n");
  DataFilter::PrintKeywordDescriptions();
}

// Exec_DataFilter::Execute()
Exec::RetType Exec_DataFilter::Execute(CpptrajState& State, ArgList& argIn) {
  DataFilter dataFilter;

  if (dataFilter.InitFilter(argIn, State.DSL(), State.DFL(), State.Debug()))
    return CpptrajState::ERR;

  size_t nelements = dataFilter.MinNumElements();
  if (nelements < 1) {
    mprinterr("Error: No data to filter. All sets must contain some data.\n");
    return CpptrajState::ERR;
  }
 
  ProgressBar progress( nelements );
  for (size_t idx = 0; idx != nelements; idx++) {
    progress.Update( idx );
    dataFilter.FilterIndex( idx );
  }

  if (dataFilter.Finalize())
    return CpptrajState::ERR;
  if (!dataFilter.IsMulti())
    mprintf("\t%u frames passed through, %u frames were filtered out.\n",
              dataFilter.Npassed(), dataFilter.Nfiltered());
  // Trigger master datafile write just in case
  State.MasterDataFileWrite();
  return CpptrajState::OK;
}

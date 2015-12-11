#include "Exec_Commands.h"
#include "CpptrajStdio.h"

void Exec_Run::Help() const {
  mprintf("  Process all trajectories currently in input trajectory list.\n"
          "  All actions in action list will be run on each frame.\n"
          "  If not processing ensemble input, all analyses in analysis\n"
          "  list will be run after trajectory processing.\n");
}
// -----------------------------------------------------------------------------
void Exec_NoExitOnError::Help() const {
  mprintf("  Do not exit when errors are encountered. This is the default\n"
          "  in interactive mode.\n");
}

Exec::RetType Exec_NoExitOnError::Execute(CpptrajState& State, ArgList&)
{
  State.SetNoExitOnError();
  mprintf("\tAttempting to ignore errors if possible.\n");
  return CpptrajState::OK;
}
// -----------------------------------------------------------------------------
void Exec_NoProgress::Help() const {
  mprintf("  Do not print progress while reading in trajectories.\n");
}

Exec::RetType Exec_NoProgress::Execute(CpptrajState& State, ArgList&)
{
  State.SetNoProgress();
  mprintf("\tProgress bar will not be used during Run.\n");
  return CpptrajState::OK;
}
// -----------------------------------------------------------------------------
void Exec_Quit::Help() const { mprintf("  Exit CPPTRAJ\n"); }
// -----------------------------------------------------------------------------
void Exec_ActiveRef::Help() const {
  mprintf("\t%s\n", DataSetList::RefArgs);
  mprintf("  Set the reference structure to be used for coordinate-based mask parsing.\n"
          "  <#> starts from 0 (first reference structure).\n");
}
// -----------------------------------------------------------------------------
void Exec_Clear::Help() const {
  mprintf("\t[ {all | <type>} ] (<type> =%s)\n", CpptrajState::PrintListKeys().c_str());
  mprintf("  Clear currently loaded objects of the specified type. If 'all' is specified\n"
          "  then clear all loaded objects.\n");
}
// -----------------------------------------------------------------------------
void Exec_RemoveData::Help() const {
  mprintf("\t[<arg>]\n"
          "  Remove data sets(s) corresponding to <arg> from data set list.\n");
}
// -----------------------------------------------------------------------------
void Exec_SetListDebug::Help() const {
  mprintf("\t[<type>] <#> (<type> =%s)\n", CpptrajState::PrintListKeys().c_str());
  mprintf("  Set debug level for new objects of the specified type. If no type is given\n"
          "  then set debug level for all new objects. Does not affect current objects.\n");
}
// -----------------------------------------------------------------------------
void Exec_ListAll::Help() const {
  mprintf("\t[<type>] (<type> =%s)\n"
          "  List currently loaded objects of the specified type. If no type is given\n"
          "  then list all loaded objects.\n", CpptrajState::PrintListKeys().c_str());
}
// -----------------------------------------------------------------------------
void Exec_SilenceActions::Help() const { mprintf("Silence Actions Init/Setup output.\n"); }
// -----------------------------------------------------------------------------
void Exec_DataFileCmd::Help() const {
  mprintf("\t<data filename> <datafile cmd>\n"
          "  Pass <datafile cmd> to specified data file currently in data file list.\n");
  DataFile::WriteHelp();
  DataFile::WriteOptions();
}


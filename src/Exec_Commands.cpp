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
// -----------------------------------------------------------------------------
void Exec_SelectAtoms::Help() const {
  mprintf("\t[%s] <mask>\n"
          "  Show atom numbers selected by <mask> for parm <parmindex>\n"
          "  (default first parm)\n", DataSetList::TopIdxArgs);
}

Exec::RetType Exec_SelectAtoms::Execute(CpptrajState& State, ArgList& argIn) {
  AtomMask tempMask( argIn.GetMaskNext() );
  Topology* parm = State.DSL()->GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  if (parm->SetupIntegerMask( tempMask )) return CpptrajState::ERR;
  mprintf("Selected %i atoms.\n", tempMask.Nselected());
  if (!argIn.hasKey("total"))
    tempMask.PrintMaskAtoms("Selected");
  return CpptrajState::OK;
}
// -----------------------------------------------------------------------------
void Exec_SelectDS::Help() const {
  mprintf("\t<dataset selection>\n"
          "  Show results of data set selection. Data set selection format is:\n"
          "\t<name>[<aspect]:<idx range>\n"
          "  Where '<name>' is the data set name, '[<aspect>]' is the data set aspect,\n"
          "  and <idx range> is a numerical range specifying data set indices (i.e. 2-5,7 etc).\n"
          "  The aspect and index portions may be optional. An asterisk '*' may be used as\n"
          "  a wildcard. E.g. 'selectds R2', 'selectds RoG[Max]', 'selectds PR[res]:2-12'\n");
}

Exec::RetType Exec_SelectDS::Execute(CpptrajState& State, ArgList& argIn) {
  std::string dsarg = argIn.GetStringNext();
  DataSetList dsets = State.DSL()->GetMultipleSets( dsarg );
  if (!dsets.empty()) {
    mprintf("SelectDS: Arg '%s':", dsarg.c_str());
    dsets.List();
  }
  return CpptrajState::OK;
}

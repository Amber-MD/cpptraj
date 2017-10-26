#include "Exec_ReadEnsembleData.h"
#include "CpptrajStdio.h"

// Exec_ReadEnsembleData::Help()
void Exec_ReadEnsembleData::Help() const
{
  mprintf("<filename> [names <additional files]\n");
}

// Exec_ReadEnsembleData::Execute()
Exec::RetType Exec_ReadEnsembleData::Execute(CpptrajState& State, ArgList& argIn)
{
  std::string additionalNames = argIn.GetStringKey("names");
  FileName fname( argIn.GetStringNext() );
  if (fname.empty()) {
    mprinterr("Error: No file name given.\n");
    return CpptrajState::ERR;
  }
  File::NameArray fileNames;
  if (!additionalNames.empty()) {
    ArgList name_list(additionalNames, ",");
    if (name_list.Nargs() < 1) {
      mprinterr("Error: No additional names specified.\n");
      return CpptrajState::ERR;
    }
    fileNames.push_back( fname );
    for (int arg = 0; arg < name_list.Nargs(); arg++)
      fileNames.push_back( name_list[arg] );
  } else {
    fileNames = File::SearchForReplicas( fname, State.Debug() );
  }
  mprintf("\t%zu files.\n", fileNames.size());
  for (File::NameArray::const_iterator it = fileNames.begin(); it != fileNames.end(); ++it)
    mprintf("\t  %s\n", it->full());
  return CpptrajState::OK;
}

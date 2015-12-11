#include "Exec_ReadData.h"
#include "CpptrajStdio.h"

void Exec_ReadData::Help() const {
  mprintf("\t<filename> [name <dsname>] [as <fmt>] [<format options>]\n"
          "  Read data from <filename> into data sets.\n");
  DataFile::ReadOptions();
}

Exec::RetType Exec_ReadData::Execute(CpptrajState& State, ArgList& argIn) {
  DataFile dataIn;
  dataIn.SetDebug( State.DFL()->Debug() );
  std::string filenameIn = argIn.GetStringNext();
  File::NameArray fnames = File::ExpandToFilenames( filenameIn );
  if (fnames.empty()) {
    mprinterr("Error: '%s' matches no files.\n", filenameIn.c_str());
    return CpptrajState::ERR;
  }
  int err = 0;
  for (File::NameArray::const_iterator fn = fnames.begin(); fn != fnames.end(); ++fn) {
    if (dataIn.ReadDataIn( *fn, argIn, *State.DSL() )!=0) {
      mprinterr("Error: Could not read data file '%s'.\n", fn->full());
      err++;
    }
  }
  if (err > 0) return CpptrajState::ERR;
  return CpptrajState::OK;
}

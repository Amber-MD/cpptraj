#include "Exec_DataFile.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString

/// Add DataSets specified by arguments to given DataFile.
// NOTE: Used by Create_DataFile and Write_DataFile
// TODO: Put in DataFile?
static int AddSetsToDataFile(DataFile& df, ArgList const& dsetArgs, DataSetList& DSL)
{
  int err = 0;
  std::string setsToWrite;
  for (ArgList::const_iterator dsa = dsetArgs.begin(); dsa != dsetArgs.end(); ++dsa) {
    DataSetList Sets = DSL.GetMultipleSets( *dsa );
    if (Sets.empty())
      mprintf("Warning: %s does not correspond to any data sets.\n", dsa->c_str());
    for (DataSetList::const_iterator set = Sets.begin(); set != Sets.end(); ++set) {
      if ( df.AddDataSet(*set) ) {
        mprinterr("Error: Could not add data set %s to file.\n", (*set)->legend());
        ++err;
      }
      setsToWrite.append( " " + (*set)->Meta().Legend() );
    }
  }
  mprintf("%s\n", setsToWrite.c_str());
  return err;
}

void Exec_CreateDataFile::Help() const {
  mprintf("\t<filename> <dataset0> [<dataset1> ...]\n"
          "  Add a file with specified data sets to the data file list. Does not\n"
          "  immediately write the data.\n");
  DataFile::WriteHelp();
  DataFile::WriteOptions();
}

Exec::RetType Exec_CreateDataFile::Execute(CpptrajState& State, ArgList& argIn)
{
  // Next string is datafile that command pertains to.
  std::string name1 = argIn.GetStringNext();
  if (name1.empty()) {
    mprinterr("Error: No filename given.\n");
    return CpptrajState::ERR;
  }
  DataFile* df = State.DFL().AddDataFile(name1, argIn);
  if (df == 0) return CpptrajState::ERR;
  return (CpptrajState::RetType)( AddSetsToDataFile(*df, argIn.RemainingArgs(), State.DSL()) );
}
// -----------------------------------------------------------------------------
void Exec_WriteDataFile::Help() const {
  mprintf("\t[<filename> <dataset0> [<dataset1> ...]]\n");
  DataFile::WriteHelp();
  mprintf("  With no arguments, write all files currently in the data file list.\n"
          "  Otherwise, write specified data sets to <filename> immediately.\n");
  DataFile::WriteOptions();
}

Exec::RetType Exec_WriteDataFile::Execute(CpptrajState& State, ArgList& argIn)
{
  // Next string is datafile that command pertains to.
  std::string name1 = argIn.GetStringNext();
  if (name1.empty()) {
    State.DFL().ResetWriteStatus();
    State.MasterDataFileWrite();
    return CpptrajState::OK;
  }
  DataFile* df = new DataFile();
  if (df == 0) return CpptrajState::ERR;
  if (!argIn.hasKey("noensextension") && State.DFL().EnsembleNum() != -1)
    name1.append( "." + integerToString(State.DFL().EnsembleNum()) );
  if (df->SetupDatafile( name1, argIn, State.Debug() )) {
    delete df;
    return CpptrajState::ERR;
  }
  mprintf("\tWriting sets to %s, format '%s'\n", df->DataFilename().full(), df->FormatString());
  int err = AddSetsToDataFile(*df, argIn.RemainingArgs(), State.DSL());
  if (err == 0) df->WriteDataOut();
  delete df;
  return (CpptrajState::RetType)err;
}

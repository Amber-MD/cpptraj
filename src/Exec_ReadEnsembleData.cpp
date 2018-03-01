#include "Exec_ReadEnsembleData.h"
#include "CpptrajStdio.h"

// Exec_ReadEnsembleData::Help()
void Exec_ReadEnsembleData::Help() const
{
  mprintf("\t<filename> [filenames <additional files>] [<readdata args>]\n"
          "  Read data sets as an ensemble, i.e. each file is a different member\n"
          "  of an ensemble. If one filename is given, it is assumed it is the\n"
          "  \"lowest\" member of an ensemble with a numerical extension, e.g.\n"
          "  file.001 and the remaining files are searched for automatically.\n"
          "  Otherwise all other members of the ensemble can be specified with\n"
          "  'filenames' and a comma-separated list e.g. \n"
          "  'file.001 filenames file.002,file.003,file.004\n"
          "  For additional 'readdata' arguments that can be passed in type\n"
          "  'help readdata'\n");
}

// Exec_ReadEnsembleData::Execute()
Exec::RetType Exec_ReadEnsembleData::Execute(CpptrajState& State, ArgList& argIn)
{
  File::NameArray fileNames;
  std::string additionalNames = argIn.GetStringKey("filenames");
  if (!additionalNames.empty()) {
    // Lowest file name specified, additional members specified with
    // filenames.
    FileName fname( argIn.GetStringNext() );
    if (fname.empty()) {
      mprinterr("Error: No file name given.\n");
      return CpptrajState::ERR;
    }
    ArgList name_list(additionalNames, ",");
    if (name_list.Nargs() < 1) {
      mprinterr("Error: No additional names specified.\n");
      return CpptrajState::ERR;
    }
    fileNames.push_back( fname );
    for (int arg = 0; arg < name_list.Nargs(); arg++)
      fileNames.push_back( name_list[arg] );
  } else {
    // Check if this corresponds to multiple file names
    fileNames = File::ExpandToFilenames( argIn.GetStringNext() );
    if (fileNames.empty()) {
      mprinterr("Error: No name or invalid filename specified.\n");
      return CpptrajState::ERR;
    }
    // If this corresponds to a single file assume this is the lowest member
    // and search for more members.
    if (fileNames.size() == 1) {
      mprintf("\tAutomatically searching for ensemble members based on filename extension.\n");
      fileNames = File::SearchForReplicas( fileNames[0], State.Debug() );
    }
  }
  mprintf("\t%zu files.\n", fileNames.size());
  for (File::NameArray::const_iterator it = fileNames.begin(); it != fileNames.end(); ++it)
    mprintf("\t  %s\n", it->full());

  unsigned int min_file = 0;
  unsigned int max_file = fileNames.size();
# ifdef MPI
  // Setup communicators if not already done.
  // NOTE: Allowing fewer threads than groups here.
  if (Parallel::SetupComms( fileNames.size(), true ))
    return CpptrajState::ERR;
  min_file = (unsigned int)Parallel::Ensemble_Beg();
  max_file = (unsigned int)Parallel::Ensemble_End();
# endif

  // Execute a data read on all files.
  int err = 0;
# ifdef MPI
  // Only thread 0 from each TrajComm does reads.
  if (!Parallel::TrajComm().Master())
    min_file = max_file + 1;
# endif
  for (unsigned int nfile = min_file; nfile < max_file; nfile++)
  {
    DataFile dataIn;
    dataIn.SetDebug( State.DFL().Debug() );
    State.DSL().SetEnsembleNum( nfile );
    // TODO ExpandToFilenames?
    if (dataIn.ReadDataIn( fileNames[nfile], argIn, State.DSL() ) != 0)
    {
      mprinterr("Error: Could not read data file '%s'\n", fileNames[nfile].full());
      err = 1;
      break;
    }
  }
# ifdef MPI
  err = Parallel::World().CheckError( err );
# endif
  if (err != 0) return CpptrajState::ERR;
# ifdef MPI
  if (Parallel::EnsembleComm().Size() > 1 && State.DFL().EnsembleNum() < 0)
    State.DFL().SetEnsembleNum( Parallel::EnsembleComm().Rank() );
# endif
  return CpptrajState::OK;
}

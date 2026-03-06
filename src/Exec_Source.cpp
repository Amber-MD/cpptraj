#include "Exec_Source.h"
#include "CpptrajStdio.h"
#include "DataIO_LeapRC.h"

// Exec_Source::Help()
void Exec_Source::Help() const
{
  mprintf("\t<filename>\n"
          " Process <filename> as a leaprc file, mostly intended for loading\n"
          " force field parameters from Amber leaprc files. If AMBERHOME\n"
          " is set CPPTRAJ will look for leaprc files in\n"
          " $AMBERHOME/dat/leap/cmd, otherwise the full path must be specified.\n");
  DataIO_LeapRC::PrintSupportedLeapCommands();
}

// Exec_Source::Execute()
Exec::RetType Exec_Source::Execute(CpptrajState& State, ArgList& argIn)
{
  // source <file>
  std::string fname = argIn.GetStringNext();
  if (fname.empty()) {
    mprinterr("Error: No filename given for 'source'.\n");
    return CpptrajState::ERR;
  }
  // If file is in this dir, read it.
  // Otherwise, prepend AMBERHOME/dat/leap/cmd
  FileName fileName;
  if ( File::Exists( fname ) ) {
    fileName.SetFileName( fname );
  } else {
    const char* env = getenv("AMBERHOME");
    if (env == 0) {
      mprinterr("Error: %s not found and AMBERHOME is not set.\n", fname.c_str());
      return CpptrajState::ERR;
    }
    fileName.SetFileName( std::string(env) + "/dat/leap/cmd/" + fname );
    if (!File::Exists(fileName)) {
      mprinterr("Error: %s not found.\n", fileName.full());
      return CpptrajState::ERR;
    }
  }
  mprintf("\tReading %s\n", fileName.full());

  DataIO_LeapRC infile;
  infile.SetDebug( State.Debug() );
  CpptrajFile tmp;
  tmp.SetupRead( fileName, State.Debug() );
  if (!infile.ID_DataFormat( tmp )) {
    mprintf("Warning: %s does not appear to be a leaprc file.\n", fileName.full());
  }
  if (infile.processReadArgs( argIn )) {
    mprinterr("Error: Could not process read args for leaprc file %s\n", fileName.full());
    return CpptrajState::ERR;
  }
  if (infile.ReadData(fileName, State.DSL(), fileName.Base())) {
    mprinterr("Error: Could not read leaprc file %s\n", fileName.full());
    return CpptrajState::ERR;
  }

  return CpptrajState::OK;
}

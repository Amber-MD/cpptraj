#include "Exec_ParseTiming.h"
#include "CpptrajStdio.h"
#include "FileName.h"
#include "BufferedLine.h"

/** CONSTRUCTOR */
Exec_ParseTiming::Exec_ParseTiming() :
  Exec(GENERAL)
{
  SetHidden( true );
}

// Exec_ParseTiming::Help()
void Exec_ParseTiming::Help() const
{

}

int Exec_ParseTiming::read_cpptraj_output(std::string const& fname) {
  BufferedLine infile;

  if (infile.OpenFileRead( fname )) {
    mprinterr("Error: Could not open '%s'\n", fname.c_str());
    return 1;
  }

  infile.CloseFile();
  return 0;
}

// Exec_ParseTiming::Execute()
Exec::RetType Exec_ParseTiming::Execute(CpptrajState& State, ArgList& argIn)
{
  File::NameArray FileNameList;

  std::string filearg = argIn.GetStringNext();
  while (!filearg.empty()) {
    File::NameArray fnames = File::ExpandToFilenames( filearg );
    if (fnames.empty()) {
      mprintf("Warning: No files matching '%s'\n", filearg.c_str());
    } else {
      for (File::NameArray::const_iterator it = fnames.begin(); it != fnames.end(); ++it)
        FileNameList.push_back( *it );
    }
    filearg = argIn.GetStringNext();
  }

  for (File::NameArray::const_iterator it = FileNameList.begin(); it != FileNameList.end(); ++it) {
    mprintf("\t%s\n", it->full());
    if (read_cpptraj_output( it->Full() )) {
      mprinterr("Error: Could not read cpptraj output '%s'\n", it->full());
      return CpptrajState::ERR;
    }
  }

  return CpptrajState::OK;
}

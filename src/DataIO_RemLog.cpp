#include "DataIO_RemLog.h"
#include "CpptrajStdio.h" 
#include "BufferedLine.h"

// CONSTRUCTOR
DataIO_RemLog::DataIO_RemLog() {}

// DataIO_RemLog::ReadData()
int DataIO_RemLog::ReadData(std::string const& fname, DataSetList& datasetlist) {
  int numexchg = 0;
  //const char* SEPARATORS = " ";
  // Buffer file
  BufferedLine buffer;
  if (buffer.OpenRead( fname )) return 1;
  buffer.SetupBuffer();

  // Read the first line. Should be '# Replica Exchange log file'
  std::string line = buffer.GetLine();
  if (line != "# Replica Exchange log file") {
    mprinterr("Error: Expected '# Replica Exchange log file', got:\n%s\n", line.c_str());
    return 1;
  }

  // Read past metadata. Save expected number of exchanges.
  while (line[0] == '#') {
    line = buffer.GetLine();
    mprintf("\t%s\n", line.c_str());
    if (line.empty()) {
      mprinterr("Error: No exchanges in rem log.\n");
      return 1;
    }
    ArgList columns( line );
    if (columns.hasKey("exchange")) break;
    if (columns.hasKey("numexchg")) {
      numexchg = columns.getNextInteger(-1);
    }
  }
  mprintf("\tRem log %s should contain %i exchanges\n", fname.c_str(), numexchg);
  if (numexchg < 1) {
    mprinterr("Error: Invalid number of exchanges (%i) in rem log.\n");
    return 1;
  }

  buffer.CloseFile();

  return 0;
}

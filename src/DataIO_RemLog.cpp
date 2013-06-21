#include "DataIO_RemLog.h"
#include "CpptrajStdio.h" 
#include "BufferedLine.h"

// CONSTRUCTOR
DataIO_RemLog::DataIO_RemLog() {}

const char* ExchgDescription[] = {
"Unknown", "Temperature", "Hamiltonian", "MultipleDim"
};

// DataIO_RemLog::ReadData()
int DataIO_RemLog::ReadData(std::string const& fname, DataSetList& datasetlist) {
  enum ExchgType { UNKNOWN = 0, TREMD, HREMD, MREMD };
  ExchgType type = UNKNOWN;
  int numexchg = 0;
  int LineNum = 0;
  //const char* SEPARATORS = " ";
  // Buffer file
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;

  // Read the first line. Should be '# Replica Exchange log file'
  std::string line = buffer.GetLine();
  ++LineNum;
  if (line.compare(0, 27, "# Replica Exchange log file") != 0) {
    mprinterr("Error: Expected '# Replica Exchange log file', got:\n%s\n", line.c_str());
    return 1;
  }

  // Read past metadata. Save expected number of exchanges.
  while (line[0] == '#') {
    line = buffer.GetLine();
    ++LineNum;
    if (line.empty()) {
      mprinterr("Error: No exchanges in rem log.\n");
      return 1;
    }
    ArgList columns( line );
    if (columns.hasKey("exchange")) break;
    mprintf("\t%s", line.c_str());
    if (columns.hasKey("numexchg")) {
      numexchg = columns.getNextInteger(-1);
    }
    if (columns.hasKey("Rep#,")) {
      if (columns[2] == "Neibr#,") type = HREMD;
      else if (columns[2] == "Velocity") type = TREMD;
    }
  }
  mprintf("\tRem log %s should contain %i exchanges\n", fname.c_str(), numexchg);
  if (numexchg < 1) {
    mprinterr("Error: Invalid number of exchanges (%i) in rem log.\n");
    return 1;
  }
  mprintf("\tExchange type is %s\n", ExchgDescription[type]);
  if (type == UNKNOWN) {
    mprinterr("Error: Could not identify exchange type.\n");
    return 1;
  }
  mprintf("\tFirst exchange is on line %i\n", LineNum);

  // Should currently be positioned at the first exchange. Need to read this
  // to determine how many replicas there are.
  int n_replicas = 0;
  const char* ptr = buffer.Line();
  while (ptr != 0 && ptr[0] != '#') {
    ptr = buffer.Line();
    ++n_replicas;
  }
  mprintf("\t%i replicas.\n", n_replicas);
  ReplicaEnsemble ensemble(n_replicas);

  // Close and reopen the file, advance back to first exchange 
  buffer.CloseFile();
  buffer.OpenFileRead( fname );
  ptr = buffer.Line();
  mprintf("FIRST LINE AFTER REOPEN: %s", ptr);
  while (ptr[0] == '#') ptr = buffer.Line();
  mprintf("%s", ptr);

  buffer.CloseFile();
  return 0;
}

#include "DataIO_CharmmRepLog.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataIO_CharmmRepLog::DataIO_CharmmRepLog()
{
  SetValid( DataSet::REMLOG );
}

// NOTE: Must match LogType
const char* DataIO_CharmmRepLog::LogDescription[] = {
  "Unknown", "Temperature", "Hamiltonian", "MultipleDim", "RXSGLD", "pH"
};

bool DataIO_CharmmRepLog::ID_DataFormat(CpptrajFile& infile) {
  // Assume file set up for read
  if (infile.OpenFile()) return false;
  // Read first two lines
  bool isrepd = false;
  ArgList line1( infile.GetLine() );
  if (line1.Nargs() == 4 &&
      line1[1] == "Replica" &&
      line1[2] == "Exchange")
  {
    const char* ptr = infile.NextLine();
    if (ptr != 0 && ptr[0] == 'R' && ptr[1] == 'E' && ptr[2] == 'X' && ptr[3] == '>')
      isrepd = true;
  }
  infile.CloseFile();
  return isrepd;
}

// DataIO_CharmmRepLog::ReadHelp()
void DataIO_CharmmRepLog::ReadHelp() {
  mprintf("\tnosearch            : Do not automatically search for MREMD dimension logs.\n"
          "\tdimfile <file>      : remd.dim file for processing MREMD logs.\n"
          "\tcrdidx <crd indices>: Use comma-separated list of indices as the initial\n"
          "\t                      coordinate indices.\n"
          "\tMultiple REM logs may be specified.\n");
}

// DataIO_CharmmRepLog::processReadArgs()
int DataIO_CharmmRepLog::processReadArgs(ArgList& argIn) {
  return 0;
}

// DataIO_CharmmRepLog::ReadData()
int DataIO_CharmmRepLog::ReadData(FileName const& fnameIn, 
                            DataSetList& datasetlist, std::string const& dsname)
{

  return 0;
}

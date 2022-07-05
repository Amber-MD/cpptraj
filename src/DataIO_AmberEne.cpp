#include "DataIO_AmberEne.h"
#include "BufferedLine.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
DataIO_AmberEne::DataIO_AmberEne()
{

}

// DataIO_AmberEne::ID_DataFormat()
bool DataIO_AmberEne::ID_DataFormat(CpptrajFile& infile)
{
  if (infile.OpenFile()) return false;
  bool isAmberEne = false;
  std::string line = infile.GetLine();
  ArgList lineArgs( line, " " );
  if (lineArgs.Nargs() > 2) {
    if (lineArgs[0] == "L0" && lineArgs[1] == "Nsteps")
      isAmberEne = true;
  }
  return isAmberEne;
}

// DataIO_AmberEne::ReadHelp()
void DataIO_AmberEne::ReadHelp()
{

}

// DataIO_AmberEne::processReadArgs()
int DataIO_AmberEne::processReadArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_AmberEne::ReadData()
int DataIO_AmberEne::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  BufferedLine infile;
  if (infile.OpenFileRead(fname)) return 1;
  const char* ptr = infile.Line();
  // Should be at first line of the header
  mprintf("DEBUG: %i '%s'\n", infile.LineNumber(), ptr);
  if (ptr[0] != 'L' || ptr[1] != '0') {
    mprinterr("Error: 'L0' not found in Amber energy file '%s'\n", fname.full());
    return 1;
  }

  typedef std::vector<std::string> Sarray;
  Sarray headers;
  unsigned int Nlines = 0;

  // Read the header
  while (ptr != 0) {
    ArgList headerArgs( ptr, " " );
    if (headerArgs.Nargs() < 1) {
      mprinterr("Error: No columns detected at line %i of Amber energy file.\n",
                infile.LineNumber());
      return 1;
    }
    if (infile.LineNumber() > 1 && headerArgs[0] == "L0") {
      // At first line of data
      break;
    }
    // Read headers
    Nlines++;
    for (int iarg = 1; iarg < headerArgs.Nargs(); iarg++)
      headers.push_back( headerArgs[iarg] );
    
    ptr = infile.Line();
  }
  mprintf("\tHeader has %zu labels over %u lines.\n", headers.size(), Nlines);
  for (Sarray::const_iterator it = headers.begin(); it != headers.end(); ++it)
    mprintf(" %s", it->c_str());
  mprintf("\n");
  return 0;
}

// DataIO_AmberEne::WriteHelp()
void DataIO_AmberEne::WriteHelp()
{

}

// DataIO_AmberEne::processWriteArgs()
int DataIO_AmberEne::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_AmberEne::WriteData()
int DataIO_AmberEne::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}

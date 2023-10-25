#include "DataIO_AmberLib.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"

/// CONSTRUCTOR
DataIO_AmberLib::DataIO_AmberLib()
{

}

// DataIO_AmberLib::ID_DataFormat()
bool DataIO_AmberLib::ID_DataFormat(CpptrajFile& infile)
{
  if (infile.OpenFile()) return false;
  std::string line = infile.GetLine();
  infile.CloseFile();
  bool isLib = (line == "!!index array str");

  return isLib;
}

// DataIO_AmberLib::ReadHelp()
void DataIO_AmberLib::ReadHelp()
{

}

// DataIO_AmberLib::processReadArgs()
int DataIO_AmberLib::processReadArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_AmberLib::ReadData()
int DataIO_AmberLib::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  BufferedLine infile;
  if (infile.OpenFileRead( fname )) {
    mprinterr("Error: Could not open Amber lib file '%s' for reading.\n", fname.full());
    return 1;
  }

  // Read first line
  std::string line = infile.GetLine();
  if (line != "!!index array str") {
    mprinterr("Error: Expected first line to be '!!index array str', got '%s'\n", line.c_str());
    return 1;
  }
  typedef std::vector<std::string> Sarray;
  Sarray UnitNames;
  // Read units
  line = infile.GetLine();
  while (!line.empty() && line[0] != '!')
  {
    UnitNames.push_back( line );
    line = infile.GetLine();
  }
  mprintf("DEBUG: Units:");
  for (Sarray::const_iterator it = UnitNames.begin(); it != UnitNames.end(); ++it)
    mprintf(" %s", it->c_str());
  mprintf("\n");

  // Now should be at first unit

  return 0;
}

// DataIO_AmberLib::WriteHelp()
void DataIO_AmberLib::WriteHelp()
{

}

// DataIO_AmberLib::processWriteArgs()
int DataIO_AmberLib::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_AmberLib::WriteData()
int DataIO_AmberLib::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}

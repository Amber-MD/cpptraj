#include "DataIO_AmberPrep.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"

/// CONSTRUCTOR
DataIO_AmberPrep::DataIO_AmberPrep()
{

}

// DataIO_AmberPrep::ID_DataFormat()
bool DataIO_AmberPrep::ID_DataFormat(CpptrajFile& infile)
{

  return false;
}

// DataIO_AmberPrep::ReadHelp()
void DataIO_AmberPrep::ReadHelp()
{

}

// DataIO_AmberPrep::processReadArgs()
int DataIO_AmberPrep::processReadArgs(ArgList& argIn)
{

  return 0;
}

static inline int CheckLine(const char* line) {
  if (line==0) {
    mprinterr("Error: Unexpected end of prep file.\n");
    return 1;
  }
  return 0;
}

// DataIO_AmberPrep::ReadData()
int DataIO_AmberPrep::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  BufferedLine infile;

  if (infile.OpenFileRead(fname)) {
    mprinterr("Error: Could not open '%s'\n", fname.full());
    return 1;
  }
  // 1 - Control for data base generation
  // IDBGEN, IREST, ITYPF
  // Format (3I)
  const char* line = infile.Line();
  if (CheckLine(line)) return 1;
  // 2 - Name of the data base file. Blank if not data base gen.
  // NAMDBF
  // Format (A80)
  line = infile.Line();
  if (CheckLine(line)) return 1;
  // 3 - Title
  // Descriptive header for the residue
  line = infile.Line();
  if (CheckLine(line)) return 1;
  mprintf("DEBUG: Prep title: '%s'\n", line);

  return 0;
}

// DataIO_AmberPrep::WriteHelp()
void DataIO_AmberPrep::WriteHelp()
{

}

// DataIO_AmberPrep::processWriteArgs()
int DataIO_AmberPrep::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_AmberPrep::WriteData()
int DataIO_AmberPrep::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}

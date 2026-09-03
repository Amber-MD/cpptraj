#include "DataIO_GBNSR6.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include <cstring>

/// CONSTRUCTOR
DataIO_GBNSR6::DataIO_GBNSR6()
{

}

// DataIO_GBNSR6::ID_DataFormat()
bool DataIO_GBNSR6::ID_DataFormat(CpptrajFile& infile)
{

  return false;
}

// DataIO_GBNSR6::ReadHelp()
void DataIO_GBNSR6::ReadHelp()
{

}

// DataIO_GBNSR6::processReadArgs()
int DataIO_GBNSR6::processReadArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_GBNSR6::ReadData()
int DataIO_GBNSR6::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  mprintf("\tReading from GBNSR6 output file %s\n", fname.base());
  BufferedLine infile;
  if (infile.OpenFileRead( fname ) != 0) {
    mprinterr("Error: Could not open file '%s'\n", fname.full());
    return 1;
  }
  const char* ptr = infile.Line();
  enum PhaseType { UNKNOWN=0, INPUT };
  PhaseType Phase = UNKNOWN;
  while (ptr != 0) {
    ArgList argline(ptr);
    if (argline.Nargs() > 0) {
      if (Phase == UNKNOWN) {
        if (argline.Nargs() >= 3 && argline[0] == "Here" && argline[1] == "is" && argline[2] == "the")
          Phase = INPUT;
      } else if (Phase == INPUT) {
        if (strncmp(ptr, "-----", 5) == 0) {
          Phase = UNKNOWN;
        } else {
          mprintf("DEBUG: [Input] %s\n", ptr);
        }
      }
    } // END nargs > 0
    ptr = infile.Line();
  }

  return 0;
}

// DataIO_GBNSR6::WriteHelp()
void DataIO_GBNSR6::WriteHelp()
{

}

// DataIO_GBNSR6::processWriteArgs()
int DataIO_GBNSR6::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_GBNSR6::WriteData()
int DataIO_GBNSR6::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}

#include "DataIO_CharmmOutput.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
DataIO_CharmmOutput::DataIO_CharmmOutput()
{

}

// DataIO_CharmmOutput::ID_DataFormat()
bool DataIO_CharmmOutput::ID_DataFormat(CpptrajFile& infile)
{
  // Assume file set up for read.
  if (infile.OpenFile()) return false;
  bool isCharmmOut = false;
  // Scan first 3 lines. That should be plenty.
  for (int num = 0; num != 3; num++)
  {
    ArgList line( infile.GetLine() );
    if (line.Nargs() >= 3 && line[0] == "Chemistry" && line[2] == "HARvard") {
      isCharmmOut = true;
      break;
    } else if (line.Nargs() >=1 && line[0] == "CHARMM>") {
      isCharmmOut = true;
      break;
    }
  }
  infile.CloseFile();

  return isCharmmOut;
}

// DataIO_CharmmOutput::ReadHelp()
void DataIO_CharmmOutput::ReadHelp()
{

}

// DataIO_CharmmOutput::WriteHelp()
void DataIO_CharmmOutput::WriteHelp()
{

}

// DataIO_CharmmOutput::processReadArgs()
int DataIO_CharmmOutput::processReadArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_CharmmOutput::ReadData()
int DataIO_CharmmOutput::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{

  return 1;
}

// DataIO_CharmmOutput::processWriteArgs()
int DataIO_CharmmOutput::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_CharmmOutput::WriteData()
int DataIO_CharmmOutput::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}

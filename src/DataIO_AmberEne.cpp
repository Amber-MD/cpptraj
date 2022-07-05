#include "DataIO_AmberEne.h"
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

  return 1;
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

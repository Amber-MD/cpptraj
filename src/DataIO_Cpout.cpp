#include <cstdio>
#include "DataIO_Cpout.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
DataIO_Cpout::DataIO_Cpout()
{

}

// DataIO_Cpout::ID_DataFormat()
bool DataIO_Cpout::ID_DataFormat(CpptrajFile& infile)
{
  bool iscpout = false;
  if (!infile.OpenFile()) {
    const char* ptr = infile.NextLine();
    if (ptr != 0) {
      float orig_ph;
      if (sscanf(ptr, "Redox potential: %f V\n", &orig_ph) == 1) {
        mprinterr("Error: Redox potential not yet supported.\n");
      } else if (sscanf(ptr, "Solvent pH: %f\n", &orig_ph) == 1) {
        ptr = infile.NextLine();
        int step_size;
        if (ptr != 0) {
          iscpout = (sscanf(ptr, "Monte Carlo step size: %d\n", &step_size) == 1);
        }
      }
    }
    infile.CloseFile();
  }
  return iscpout;
}

// DataIO_Cpout::ReadHelp()
void DataIO_Cpout::ReadHelp()
{

}

// DataIO_Cpout::processReadArgs()
int DataIO_Cpout::processReadArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_Cpout::ReadData()
int DataIO_Cpout::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{

  return 1;
}

// DataIO_Cpout::WriteHelp()
void DataIO_Cpout::WriteHelp()
{

}

// DataIO_Cpout::processWriteArgs()
int DataIO_Cpout::processWriteArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_Cpout::WriteData()
int DataIO_Cpout::WriteData(FileName const& fname, DataSetList const& dsl)
{

  return 1;
}

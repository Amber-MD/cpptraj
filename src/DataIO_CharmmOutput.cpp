#include "DataIO_CharmmOutput.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"

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

// DataIO_CharmmOutput::processReadArgs()
int DataIO_CharmmOutput::processReadArgs(ArgList& argIn)
{

  return 0;
}

// DataIO_CharmmOutput::ReadData()
int DataIO_CharmmOutput::ReadData(FileName const& fname, DataSetList& dsl, std::string const& dsname)
{
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;
  const char* ptr = buffer.Line();
  // Need to scan down until we get to DYNA
  bool DYNA = false;
  while (ptr != 0) {
    if ( ptr[0] == 'D' && ptr[1] == 'Y' && ptr[2] == 'N' && ptr[3] == 'A' ) {
      DYNA = true;
      break;
    }
    ptr = buffer.Line();
  }
  if (!DYNA) {
    mprinterr("Error: 'DYNA' not found in output file '%s'.\n", fname.full());
    return 1;
  }
  // Figure out what terms we have. Format is 'DYNA <Type>: E0 E1 ...'
  typedef std::vector<std::string> Sarray;
  Sarray Terms;
  while (ptr != 0 && ptr[1] != '-') {
    ArgList line( ptr );
    for (int col = 2; col < line.Nargs(); col++)
      Terms.push_back( line[col] );
    ptr = buffer.Line();
  }
  mprintf("\t%zu terms:", Terms.size());
  for (Sarray::const_iterator it = Terms.begin(); it != Terms.end(); ++it)
    mprintf(" %s", it->c_str());
  mprintf("\n");
  
  return 0;
}

// DataIO_CharmmOutput::WriteHelp()
void DataIO_CharmmOutput::WriteHelp()
{

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

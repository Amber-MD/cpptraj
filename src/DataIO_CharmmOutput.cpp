#include <cstdlib> // atof
#include "DataIO_CharmmOutput.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "DataSet_double.h"

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
  // CHARMM is inconsistent in that the 'DYNA DYN:' header line does not
  // match the corresponding 'DYNA>' line; it is missing 'DYN>'. This
  // is OK since we do not really care about the Step anyway, but
  // this means 'Step' is skipped so do not create a DataSet for it.
  typedef std::vector<std::string> Sarray;
  Sarray Terms;
  int timeIdx = -1;
  while (ptr != 0 && ptr[1] != '-') {
    ArgList line( ptr );
    for (int col = 2; col < line.Nargs(); col++) {
      if (line[col] == "Time") timeIdx = (int)Terms.size();
      if (line[col] != "Step")
        Terms.push_back( line[col] );
    }
    ptr = buffer.Line();
  }
  mprintf("\t%zu terms:", Terms.size());
  DataSetList::DataListType inputSets;
  inputSets.reserve( Terms.size() );
  for (Sarray::const_iterator it = Terms.begin(); it != Terms.end(); ++it) {
    mprintf(" %s", it->c_str());
    inputSets.push_back( new DataSet_double() );
    inputSets.back()->SetMeta( MetaData(dsname, it->substr(0, 4)) );
  }
  mprintf("\n");
  mprintf("DEBUG: Time index: %i\n", timeIdx);
  mprintf("DEBUG: [%s]\n", ptr);

  // Read data
  int step = 0;
  bool readFile = true;
  while (readFile) {
    // Scan to next DYNA> section
    while (ptr != 0 &&
           ptr[0] != 'D' && ptr[1] != 'Y' && ptr[2] != 'N' && ptr[3] != 'A' && ptr[4] != '>' )
      ptr = buffer.Line();
    mprintf("%i [%s]\n", buffer.LineNumber(), ptr);
    if (ptr == 0)
      readFile = false;
    else {
      int idx = 0;
      while (ptr != 0 && ptr[1] != '-') {
        int ntoken = buffer.TokenizeLine(" ");
        buffer.NextToken(); // DYNA
        buffer.NextToken(); // <type>
        for (int token = 2; token < ntoken; token++) {
          double dval = atof( buffer.NextToken() );
          mprintf("DEBUG: [%i] '%s' %12.5f\n", idx, Terms[idx].c_str(), dval);
          inputSets[idx++]->Add(step, &dval);
        }
        ptr = buffer.Line();
      }
    }
  }
  mprintf("\n");
  buffer.CloseFile();

  
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

#include <cstdio> // sscanf
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
  // Terms start at character 14.
  typedef std::vector<std::string> Sarray;
  Sarray Terms;
  Sarray LineHeaders;
  int timeIdx = -1;
  typedef std::vector<int> Iarray;
  Iarray nTermsInLine;
  // CHARMM data does not appear to have consistent formatting,
  // e.g. in the XTLE> line numbers can start at the 2nd data column, so
  // need to determine start column.
  Iarray startColumn;
  while (ptr != 0 && ptr[1] != '-') {
    // Determine line header. DYN gets no header.
    std::string headerLine(ptr+5, 8);
    ArgList header(headerLine, " :");
    mprintf("DEBUG: header '%s'\n", header[0].c_str());
    ArgList line( ptr+14 );
    for (int col = 0; col < line.Nargs(); col++) {
      if (line[col] == "Time") timeIdx = (int)Terms.size();
    }
    size_t c0 = line.ArgLineStr().find( line[0] );
    startColumn.push_back( (int)(c0 / 13) );
    mprintf("DEBUG: c0= %zu\n", c0);
    // Next line
    ptr = buffer.Line();
    nTermsInLine.push_back( line.Nargs() );
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
  mprintf("DEBUG: %zu lines:\n", nTermsInLine.size());
  for (unsigned int i = 0; i != nTermsInLine.size(); i++)
    mprintf("DEBUG:\t\t[%u] %i terms, start col %i.\n", i+1, nTermsInLine[i], startColumn[i]);
/*
  // Read data.
  int step = 0;
  bool readFile = true;
  while (readFile) {
    // Scan to next DYNA> section
    while (ptr != 0) {
      if (ptr[0] == 'D' && ptr[1] == 'Y' && ptr[2] == 'N' && ptr[3] == 'A' && ptr[4] == '>' )
        break;
      ptr = buffer.Line();
    }
    if (ptr == 0)
      readFile = false;
    else {
      double dvals[5];
      std::fill(dvals, dvals+5, 0.0); // DEBUG
      mprintf("%i [%s]\n", buffer.LineNumber(), ptr);
      for (unsigned int i = 0; i != nTermsInLine.size(); i++) {
        int nvals = sscanf(ptr+14,"%13lf%13lf%13lf%13lf%13lf",
                           dvals, dvals+1, dvals+2, dvals+3, dvals+4);
        // SANITY CHECK
        if (nvals != nTermsInLine[i]) {
          mprinterr("Error: Number of terms in line %i (%i) != expected terms (%i)\n",
                    buffer.LineNumber(), nvals, nTermsInLine[i]);
          return 1; // TODO carry on?
        }
        mprintf("DEBUG:");
        for (int val = 0; val < nvals; val++)
          mprintf(" %13.5f", dvals[val]);
        mprintf("\n");
        ptr = buffer.Line();
      }
    }
  }
*/
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

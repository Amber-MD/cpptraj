#include <cstdio> // sscanf
#include "DataIO_CharmmOutput.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "DataSet_double.h"

/// CONSTRUCTOR
DataIO_CharmmOutput::DataIO_CharmmOutput() { }

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
void DataIO_CharmmOutput::ReadHelp() { }

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
  while (ptr != 0 && ptr[1] != '-') {
    // Determine line header.
    // NOTE DYN gets no header in DYNA> section.
    std::string headerLine(ptr+5, 8);
    ArgList header(headerLine, " :");
    LineHeaders.push_back( header[0] + ">" );
    ArgList line( ptr+14 );
    for (int col = 0; col < line.Nargs(); col++) {
      if (line[col] == "Time") timeIdx = (int)Terms.size();
      Terms.push_back( line[col] );
    }
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
    inputSets.back()->SetMeta( MetaData(dsname, *it) );
  }
  mprintf("\n");
  if (debug_ > 0) {
    mprintf("DEBUG: Time index: %i\n", timeIdx);
    mprintf("DEBUG: [%s]\n", ptr);
    mprintf("DEBUG: %zu lines:\n", nTermsInLine.size());
    for (unsigned int i = 0; i != nTermsInLine.size(); i++)
      mprintf("DEBUG:\t\t[%u] '%s' %i terms.\n",
              i+1, LineHeaders[i].c_str(), nTermsInLine[i]);
  }
  // Read data.
  int set = 0;
  int lastStep = -1;
  bool readFile = true;
  bool isRepD = false;
  while (readFile) {
    // Scan to next DYNA> section
    while (ptr != 0) {
      if (ptr[0] == 'D' && ptr[1] == 'Y' && ptr[2] == 'N' && ptr[3] == 'A' && ptr[4] == '>' )
        break;
      if (ptr[0] == 'R' && ptr[1] == 'E' && ptr[2] == 'P' && ptr[3] == 'D' && ptr[4] == '>' )
        isRepD = true;
      ptr = buffer.Line();
    }
    if (ptr == 0)
      readFile = false;
    else {
      int step;
      double dvals[5];
      bool ignoreStep = false;
      std::fill(dvals, dvals+5, 0.0); // DEBUG
      //mprintf("%i [%s]\n", buffer.LineNumber(), ptr);
      int idx = 0; // Index into inputSets
      for (unsigned int i = 0; i != nTermsInLine.size(); i++) {
        // Determine if this line is present. First line should always
        // be present, DYN
        bool lineIsPresent = true;
        if (i > 0) {
          lineIsPresent = (LineHeaders[i].compare(0, LineHeaders[i].size(),
                                                  ptr+5, LineHeaders[i].size()) == 0);
        } else {
          // Line 0 should have the step
          sscanf(ptr+5, "%9i", &step);
          //mprintf("DEBUG: Step %i\n", step);
          if (step == lastStep) {
            // If REPD, dynamics are restarted after exchange, so the output
            // from the last step is repeated and should be ignored.
            if (isRepD)
              ignoreStep = true;
            else
              mprintf("Warning: Repeated step detected in non-REPD run: %i\n", step);
          }
        }
        if (ignoreStep) {
          ptr = buffer.Line();
        } else if (lineIsPresent) {
          int nvals = sscanf(ptr+14,"%13lf%13lf%13lf%13lf%13lf",
                             dvals, dvals+1, dvals+2, dvals+3, dvals+4);
          // SANITY CHECK
          if (nvals != nTermsInLine[i]) {
            mprinterr("Error: Number of terms in line %i (%i) != expected terms (%i)\n",
                      buffer.LineNumber(), nvals, nTermsInLine[i]);
            return 1; // TODO carry on?
          }
          //mprintf("DEBUG: %8s", LineHeaders[i].c_str());
          for (int val = 0; val < nvals; val++) {
            //mprintf(" %13.5f", dvals[val]);
            inputSets[idx+val]->Add(set, dvals + val);
          }
          //mprintf("\n");
          ptr = buffer.Line();
        }
        idx += nTermsInLine[i];
      } // END loop over DYNA> lines
      if (step != lastStep) set++;
      lastStep = step;
    } // END start DYNA> section
  }
  if (isRepD)
    mprintf("\tREPD run detected. Repeated steps were ignored.\n");
  buffer.CloseFile();
  // Separate out time values.
  
  DataSetList::DataListType dataSets;
  dataSets.reserve( inputSets.size() - 1 );
  for (int idx = 0; idx != (int)inputSets.size(); idx++) {
    if (idx != timeIdx)
      dataSets.push_back( inputSets[idx] );
  }
  DataSetList::Darray const& timeVals = ((DataSet_double*)inputSets[timeIdx])->Data();
  if (dsl.AddOrAppendSets( "Time", timeVals, dataSets )) return 1;
  
  return 0;
}

// DataIO_CharmmOutput::WriteHelp()
void DataIO_CharmmOutput::WriteHelp() { }

// DataIO_CharmmOutput::processWriteArgs()
int DataIO_CharmmOutput::processWriteArgs(ArgList& argIn)
{
  return 0;
}

// DataIO_CharmmOutput::WriteData()
int DataIO_CharmmOutput::WriteData(FileName const& fname, DataSetList const& dsl)
{
  mprinterr("Error: Charmm output write not supported.\n");
  return 1;
}

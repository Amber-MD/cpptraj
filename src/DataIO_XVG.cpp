#include <cstdlib> // atof
#include "DataIO_XVG.h"
#include "BufferedLine.h"
#include "CpptrajStdio.h"
#include "DataSet_double.h"

bool DataIO_XVG::ID_DataFormat(CpptrajFile& infile) {
  if (infile.OpenFile()) return false;
  const char* ptr = infile.NextLine();
  while (ptr != 0 && ptr[0] == '#') {
    const char* cc = ptr;
    while (*cc != '\0') {
      if (*cc == 'G') {
        if ( cc[2] == 'R' && cc[4]  == 'O' && cc[6]  == 'M' &&
             cc[8] == 'A' && cc[10] == 'A' && cc[12] == 'C' )
        {
          infile.CloseFile();
          mprintf("DEBUG:\tFound G R O M A C\n");
          return true;
        }
      }
      ++cc;
    }
    ptr = infile.NextLine();
  }
  infile.CloseFile();
  return false;
}

int DataIO_XVG::ReadData(FileName const& fname, 
                         DataSetList& datasetlist, std::string const& dsname)
{
  std::vector<std::string> Legends;
  BufferedLine infile;

  if (infile.OpenFileRead( fname )) return 1;
  const char* ptr = infile.Line();
  if (ptr == 0) return 1;
  // Skip any comments
  while (ptr != 0 && ptr[0] == '#')
    ptr = infile.Line();
  // Try to get set legends
  while (ptr != 0 && ptr[0] == '@') {
    ArgList line(ptr, " \t");
    if (line.Nargs() > 3 && line[1][0] == 's') {
      std::string legend = line.GetStringKey("legend");
      if (!legend.empty()) {
        // Spaces will cause issues with data set selection.
        for (std::string::iterator s = legend.begin(); s != legend.end(); ++s)
          if (*s == ' ') *s = '_';
        Legends.push_back( legend );
      }
    }
    ptr = infile.Line();
  }
  if (Legends.empty()) {
    mprintf("Warning: No set legends found in XVG file, assuming only one data set.\n");
  }
  if (ptr == 0) {
    mprinterr("Error: No data in XVG file.\n");
    return 1;
  }
  // Create 1 data set for each legend
  DataSetList::DataListType inputSets;
  if (!Legends.empty()) {
    for (unsigned int i = 0; i != Legends.size(); i++) {
      MetaData md( dsname, i );
      md.SetLegend( Legends[i] );
      DataSet_double* ds = new DataSet_double();
      if (ds == 0) return 1;
      ds->SetMeta( md );
      inputSets.push_back( ds );
    }
  } else {
    MetaData md(dsname);
    DataSet_double* ds = new DataSet_double();
    if (ds == 0) return 1;
    ds->SetMeta( md );
    inputSets.push_back( ds );
  }
  mprintf("\t%s has %zu columns of data.\n", fname.base(), inputSets.size());
  // Should now be positioned at first line of data. Assume first column is time values.
  DataSetList::Darray Xvals;
  int expectedCols = (int)inputSets.size() + 1;
  while (ptr != 0) {
    int ncols = infile.TokenizeLine(" \t");
    if (ncols != expectedCols)
      mprinterr("Error: Line %i: %i columns != expected # cols %i\n", infile.LineNumber(),
                ncols, expectedCols);
    else {
      Xvals.push_back( atof( infile.NextToken() ) );
      for (unsigned int i = 0; i != inputSets.size(); i++)
        ((DataSet_double*)inputSets[i])->AddElement( atof( infile.NextToken() ) );
    }
    ptr = infile.Line();
  }
  infile.CloseFile();
  return (datasetlist.AddOrAppendSets( "", Xvals, inputSets ));
}

#include <cstdlib> // atof
#include "DataIO_XVG.h"
#include "BufferedLine.h"
#include "CpptrajStdio.h"

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

int DataIO_XVG::ReadData(std::string const& fname, 
                         DataSetList& datasetlist, std::string const& dsname)
{
  typedef std::vector<double> Darray;
  std::vector<Darray> Dsets;
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
    mprinterr("Error: No set legends found in XVG file.\n");
    return 1;
  }
  if (ptr == 0) {
    mprinterr("Error: No data in XVG file.\n");
    return 1;
  }
  Dsets.resize( Legends.size() + 1); // +1 for x values
  mprintf("\t%s has %zu columns of data (including time values).\n", infile.Filename().base(),
          Dsets.size());
  // Should now be positioned at first line of data
  while (ptr != 0) {
    int ncols = infile.TokenizeLine(" \t");
    if (ncols != (int)Dsets.size())
      mprinterr("Error: Line %i: %i columns != expected # sets %zu\n", infile.LineNumber(),
                ncols, Dsets.size());
    else {
      for (std::vector<Darray>::iterator set = Dsets.begin(); set != Dsets.end(); ++set)
        set->push_back( atof( infile.NextToken() ) );
    }
    ptr = infile.Line();
  }
  // Save data sets
  for (unsigned int i = 1; i != Dsets.size(); i++) {
    DataSet* ds = datasetlist.AddOrAppendSet(dsname, i-1, Legends[i-1], Dsets[0], Dsets[i]);
    if (ds == 0) return 1;
    ds->SetLegend( Legends[i-1] );
  }

  infile.CloseFile();
  return 0;
}

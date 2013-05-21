#include <cstdlib> // atoi, atof
#include <cstring> // strchr
#include <cctype>  // isdigit, isalpha
#include "DataIO_Std.h"
#include "CpptrajStdio.h" 
#include "StringRoutines.h" // SetStringFormatString
#include "BufferedLine.h"
#include "Array1D.h"
#include "DataSet_2D.h"

// CONSTRUCTOR
DataIO_Std::DataIO_Std() : 
  hasXcolumn_(true), 
  writeHeader_(true), 
  square2d_(false) {}

static void PrintColumnError() {
  mprinterr("Error: Number of columns in file changes.\n");
}

// DataIO_Std::ReadData()
int DataIO_Std::ReadData(std::string const& fname, DataSetList& datasetlist) {
  ArgList labels;
  bool hasLabels = false;
  Array1D DsetList;
  int indexcol = -1;
  const char* SEPARATORS = " ,\t"; // whitespace, comma, or tab-delimited
  //DataSet::DataType indextype = DataSet::UNKNOWN_DATA;

  // Buffer file
  BufferedLine buffer;
  if (buffer.OpenRead( fname )) return 1;
  buffer.SetupBuffer();

  // Read the first line. Attempt to determine the number of columns
  const char* linebuffer = buffer.Line();
  if (linebuffer == 0) return 1;
  int ntoken = buffer.TokenizeLine( SEPARATORS );
  if ( ntoken == 0 ) {
    mprinterr("Error: No columns detected in %s\n", buffer.Filename().full());
    return 1;
  }

  // If first line begins with a '#', assume it contains labels
  if (linebuffer[0]=='#') {
    labels.SetList(linebuffer+1, SEPARATORS );
    hasLabels = true;
    // If label is Frame assume it is the index column
    if (labels[0] == "Frame") 
      indexcol = 0;
    // Read in next non # line, should be data.
    while (linebuffer[0] == '#') {
      linebuffer = buffer.Line();
      if (linebuffer == 0) return 1;
    }
    if (buffer.TokenizeLine( SEPARATORS ) != ntoken) {
      PrintColumnError();
      return 1;
    }
  }

  // Determine the type of data stored in each column 
  for (int col = 0; col < ntoken; ++col) {
    const char* token = buffer.NextToken();
    //mprintf("\t\tFirst Data in col %i = %s", col+1, token);
    //if (col == indexcol)
    //  mprintf(" INDEX");
    // Determine data type
    DataSet_1D* dset = 0;
    if ( isdigit( token[0] )    || 
                  token[0]=='+' || 
                  token[0]=='-' ||
                  token[0]=='.'   )
    {
      if ( strchr( token, '.' ) != 0 ) {
        //mprintf(" DOUBLE!\n");
        if ( col != indexcol )
          dset = (DataSet_1D*)datasetlist.AddSetIdx( DataSet::DOUBLE, buffer.Filename().Base(), col+1 );
        //else
        //  indextype = DataSet::DOUBLE;
      } else {
        //mprintf(" INTEGER!\n");
        if (col != indexcol)
          dset = (DataSet_1D*)datasetlist.AddSetIdx( DataSet::INTEGER, buffer.Filename().Base(), col+1 );
        //else
        //  indextype = DataSet::INT;
      }
    } else {
      // Assume string
      //mprintf(" STRING!\n");
      // STRING columns cannot be index columns
      if ( col == indexcol ) {
        mprinterr("Error: DataFile %s index column %i has string values.\n", 
                  buffer.Filename().full(), indexcol+1);
        return 1;
      }
      dset = (DataSet_1D*)datasetlist.AddSetIdx( DataSet::STRING, buffer.Filename().Base(), col+1 );
    } 
    // Set legend to label if present
    if ( dset != 0 && hasLabels)
      dset->SetLegend( labels[col] );

    if ( col != indexcol && dset == 0 ) {
      mprinterr("Error: DataFile %s: Could not identify column %i", 
                buffer.Filename().full(), col+1);
      mprinterr(" (token=%s)\n",token);
      return 1;
    }
    DsetList.push_back( dset );
  }

  int ival = 0;
  double dval = 0;
  // NOTE: Temporarily disabling index.
  int indexval = -1; // So that when no index col first val incremented to 0
  do {
    if ( buffer.TokenizeLine( SEPARATORS ) != ntoken ) {
      PrintColumnError();
      break;
    } 
    // Deal with index.
/*    if (indexcol != -1) {
      switch ( indextype ) {
        case DataSet::INT   : indexval = convertToInteger( dataline[indexcol] ); break;
        case DataSet::DOUBLE: indexval = (int)convertToDouble( dataline[indexcol] ); break;
        default: return 1;
      }
      // FIXME: Subtracting 1 since everything should go from 0
      --indexval;
    } else {*/
      ++indexval;
//    }
    // Convert data in columns
    for (int i = 0; i < ntoken; ++i) {
      const char* token = buffer.NextToken();
      if (DsetList[i] == 0) continue;
      switch ( DsetList[i]->Type() ) {
        case DataSet::INTEGER: 
          ival = atoi( token ); 
          DsetList[i]->Add( indexval, &ival );
          break;
        case DataSet::DOUBLE: 
          dval = atof( token ); 
          DsetList[i]->Add( indexval, &dval );
          break;
        case DataSet::STRING: 
          DsetList[i]->Add( indexval, (char*)token ); // TODO: Fix cast
          break;
        default: continue; 
      }
    }
  } while (buffer.Line() != 0);
  buffer.CloseFile();
  mprintf("\tDataFile %s has %i columns.\n", buffer.Filename().full(), ntoken);
  if (hasLabels) {
    mprintf("\tDataFile contains labels:\n");
    labels.PrintList();
  }
  if (indexcol != -1)
    mprintf("\tIndex column is %i\n", indexcol + 1);

  return 0;
}

// DataIO_Std::processWriteArgs()
int DataIO_Std::processWriteArgs(ArgList &argIn) {
  hasXcolumn_ = !argIn.hasKey("noxcol");
  writeHeader_ = !argIn.hasKey("noheader");
  square2d_ = argIn.hasKey("square2d");
  return 0;
}

// -----------------------------------------------------------------------------
// WriteNameToBuffer()
static void WriteNameToBuffer(CpptrajFile& fileIn, DataSet const& DS, bool leftAlign) 
{
  std::string temp_name = DS.Legend();
  // If left aligning, add '#' to name; ensure that name will not be
  // larger than column width.
  if (leftAlign) {
    if (temp_name[0]!='#')
      temp_name.insert(0,"#");
  }
  int width = DS.ColumnWidth();
  if ((int)temp_name.size() > width)
    temp_name.resize( width );
  // Replace any spaces with underscores
  for (std::string::iterator tc = temp_name.begin(); tc != temp_name.end(); ++tc)
    if ( *tc == ' ' )
      *tc = '_';
  // Set up header format string
  std::string header_format = SetStringFormatString(width, leftAlign);
  fileIn.Printf(header_format.c_str(), temp_name.c_str());
}
// -----------------------------------------------------------------------------

// DataIO_Std::WriteData()
int DataIO_Std::WriteData(std::string const& fname, DataSetList const& SetList,
                          DimArray const& Dim) 
{
  std::string x_col_format, x_label;
  int xcol_width = 8;
  int xcol_precision = 3;

  // Hold all 1D data sets.
  Array1D Sets( SetList );
  if (Sets.empty()) return 1;
  // Use X dimension of set 0 for all set dimensions.
  // TODO: Check for empty dim.
  Dimension const& Xdim = Dim[0];

  // Determine size of largest DataSet.
  size_t maxFrames = Sets.DetermineMax();

  // Set up X column.
  if (hasXcolumn_) {
    // Set up x column label. Precede label with a '#'.
    x_label = Xdim.Label();
    x_label.insert(0, "#");
    if ((int)x_label.size() > xcol_width) xcol_width = (int)x_label.size();
    // Create format string for X column based on dimension in first data set.
    if (Xdim.Step() == 1.0) xcol_precision = 0;
    x_col_format = SetupCoordFormat( maxFrames, Xdim, xcol_width, xcol_precision ); 
  } else {
    // If not writing an X-column, set the format for the first dataset
    // to left-aligned.
    if (!hasXcolumn_) Sets[0]->SetDataSetFormat( true );
  }

  // Open output file.
  CpptrajFile file;
  if (file.OpenWrite( fname )) return 1;

  // Write header to buffer
  if (writeHeader_) {
    // If x-column present, write x-label
    if (hasXcolumn_) {
      std::string x_header_fmt = SetStringFormatString( xcol_width, true );
      file.Printf(x_header_fmt.c_str(), x_label.c_str());
    }
    // Write dataset names to header, left-aligning first set if no X-column
    Array1D::const_iterator set = Sets.begin();
    if (!hasXcolumn_)
      WriteNameToBuffer( file, *(*set), true );
    else
      WriteNameToBuffer( file, *(*set), false );
    ++set;
    for (; set != Sets.end(); ++set) 
      WriteNameToBuffer( file, *(*set), false );
    file.Printf("\n"); 
  }

  // Write Data
  for (size_t frame=0L; frame < maxFrames; frame++) {
    // Output Frame for each set
    if (hasXcolumn_)
      file.Printf(x_col_format.c_str(), Xdim.Coord(frame));
    for (Array1D::const_iterator set = Sets.begin(); set != Sets.end(); ++set) 
      (*set)->WriteBuffer(file, frame);
    file.Printf("\n"); 
  }
  file.CloseFile();
  return 0;
}

// DataIO_Std::WriteDataInverted()
int DataIO_Std::WriteDataInverted(std::string const& fname, DataSetList const& SetList,
                                  DimArray const& Dim)
{
  // Hold all 1D data sets.
  Array1D Sets( SetList );
  if (Sets.empty()) return 1;
  // Determine size of largest DataSet.
  size_t maxFrames = Sets.DetermineMax();
  // Open output file
  CpptrajFile file;
  if (file.OpenWrite( fname )) return 1;
  // Write each set to a line.
  for (Array1D::const_iterator set = Sets.begin(); set != Sets.end(); ++set) {
    // Write dataset name as first column.
    WriteNameToBuffer( file, *(*set), false); 
    // Write each frame to subsequent columns
    for (size_t frame=0L; frame < maxFrames; frame++) 
      (*set)->WriteBuffer(file, frame);
    file.Printf("\n");
  }
  file.CloseFile();
  return 0;
}

// DataIO_Std::WriteData2D()
int DataIO_Std::WriteData2D( std::string const& fname, DataSet const& setIn, 
                             DimArray const& Dim)
{
  if (setIn.Ndim() != 2) {
    mprinterr("Internal Error: DataSet %s in DataFile %s has %zu dimensions, expected 2.\n",
              setIn.Legend().c_str(), fname.c_str(), setIn.Ndim());
    return 1;
  }
  DataSet_2D const& set = static_cast<DataSet_2D const&>( setIn );
  // Open output file
  CpptrajFile file;
  if (file.OpenWrite( fname )) return 1;
  
  if (square2d_) {
    std::string ycoord_fmt;
    // Print XY values in a grid:
    // x0y0 x1y0 x2y0
    // x0y1 x1y1 x2y1
    // x0y2 x1y2 x2y2
    // If file has header, top-left value will be '#<Xlabel>-<Ylabel>',
    // followed by X coordinate values.
    if (writeHeader_) {
      std::string header = "#" + Dim[0].Label() + "-" + Dim[1].Label();
      int col1_width = (int)header.size();
      std::string col1_fmt = SetStringFormatString( col1_width, true );
      file.Printf(col1_fmt.c_str(), header.c_str());
      std::string xcoord_fmt = SetupCoordFormat( set.Ncols(), Dim[0], set.ColumnWidth(), 3 );
      for (size_t ix = 0; ix < set.Ncols(); ix++)
        file.Printf(xcoord_fmt.c_str(), Dim[0].Coord( ix ));
      file.Printf("\n");
      ycoord_fmt = SetupCoordFormat( set.Nrows(), Dim[1], col1_width, 3 );
    }
    for (size_t iy = 0; iy < set.Nrows(); iy++) {
      if (writeHeader_)
        file.Printf(ycoord_fmt.c_str(), Dim[1].Coord( iy ));
      for (size_t ix = 0; ix < set.Ncols(); ix++)
        set.Write2D( file, ix, iy );
      file.Printf("\n");
    }
  } else {
    // Print X Y Values
    // x y val(x,y)
    if (writeHeader_)
      file.Printf("#%s %s %s\n", Dim[0].Label().c_str(), 
                  Dim[1].Label().c_str(), set.Legend().c_str());
    std::string col_fmt = SetupCoordFormat( set.Ncols(), Dim[0], 8, 3 ) + " " +
                          SetupCoordFormat( set.Nrows(), Dim[1], 8, 3 ); 
    for (size_t iy = 0; iy < set.Nrows(); ++iy) {
      for (size_t ix = 0; ix < set.Ncols(); ++ix) {
        file.Printf(col_fmt.c_str(), Dim[0].Coord( ix ), Dim[1].Coord( iy ));
        set.Write2D( file, ix, iy );
        file.Printf("\n");
      }
    }
  }
  return 0;
}

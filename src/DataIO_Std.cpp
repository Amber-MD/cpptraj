#include <cstdio>  // sscanf
#include <cstdlib> // atoi, atof
#include <cstring> // strchr
#include <cctype>  // isdigit, isalpha
#include <cmath>   // modf TODO put function in StringRoutines?
#include "DataIO_Std.h"
#include "CpptrajStdio.h" 
#include "StringRoutines.h" // SetStringFormatString
#include "BufferedLine.h"
#include "TextFormat.h"
#include "DataSet_double.h" // For reading TODO remove dependency?
#include "DataSet_string.h" // For reading TODO remove dependency?
#include "DataSet_Vector.h" // For reading TODO remove dependency?
#include "DataSet_Mat3x3.h" // For reading TODO remove dependency?
#include "DataSet_2D.h"
#include "DataSet_3D.h"
#include "DataSet_Cmatrix.h"

// CONSTRUCTOR
DataIO_Std::DataIO_Std() :
  DataIO(true, true, true), // Valid for 1D, 2D, 3D
  mode_(READ1D),
  indexcol_(-1),
  isInverted_(false), 
  hasXcolumn_(true), 
  writeHeader_(true), 
  square2d_(false)
{}

static void PrintColumnError(int idx) {
  mprinterr("Error: Number of columns in file changes at line %i.\n", idx);
}

void DataIO_Std::ReadHelp() {
  mprintf("\tread1d:      Read data as 1D data sets (default).\n"
          "\tread2d:      Read data as 2D square matrix.\n"
          "\tvector:      Read data as vector: VX VY VZ [OX OY OZ]\n"
          "\tmat3x3:      Read data as 3x3 matrices: M(1,1) M(1,2) ... M(3,2) M(3,3)\n"
          "\tindex <col>: (1D) Use column # (starting from 1) as index (X) column.\n");

}

const char* DataIO_Std::SEPARATORS = " ,\t"; // whitespace, comma, or tab-delimited

int DataIO_Std::processReadArgs(ArgList& argIn) {
  mode_ = READ1D;
  if (argIn.hasKey("read1d")) mode_ = READ1D;
  else if (argIn.hasKey("read2d")) mode_ = READ2D;
  else if (argIn.hasKey("vector")) mode_ = READVEC;
  else if (argIn.hasKey("mat3x3")) mode_ = READMAT3X3;
  indexcol_ = argIn.getKeyInt("index", -1);
  // Column user args start from 1.
  if (indexcol_ == 0) {
    mprinterr("Error: Column numbering for standard data files starts from 1.\n");
    return 1;
  }
  if (indexcol_ > 0) --indexcol_;
  return 0;
}
  
// TODO: Set dimension labels
// DataIO_Std::ReadData()
int DataIO_Std::ReadData(FileName const& fname, 
                         DataSetList& dsl, std::string const& dsname)
{
  int err = 0;
  switch ( mode_ ) {
    case READ1D: err = Read_1D(fname.Full(), dsl, dsname); break;
    case READ2D: err = Read_2D(fname.Full(), dsl, dsname); break;
    case READVEC: err = Read_Vector(fname.Full(), dsl, dsname); break;
    case READMAT3X3: err = Read_Mat3x3(fname.Full(), dsl, dsname); break;
  }
  return err;
}

// DataIO_Std::Read_1D()
int DataIO_Std::Read_1D(std::string const& fname, 
                        DataSetList& datasetlist, std::string const& dsname)
{
  ArgList labels;
  bool hasLabels = false;
  // Buffer file
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;

  // Read the first line. Attempt to determine the number of columns
  const char* linebuffer = buffer.Line();
  if (linebuffer == 0) return 1;
  int ntoken = buffer.TokenizeLine( SEPARATORS );
  if ( ntoken == 0 ) {
    mprinterr("Error: No columns detected in %s\n", buffer.Filename().full());
    return 1;
  }

  // Try to skip past any comments. If line begins with a '#', assume it
  // contains labels. 
  bool isCommentLine = true;
  const char* ptr = linebuffer;
  while (isCommentLine) {
    // Skip past any whitespace
    while ( *ptr != '\0' && isspace(*ptr) ) ++ptr;
    // Assume these are column labels until proven otherwise.
    if (*ptr == '#') {
      labels.SetList(ptr+1, SEPARATORS );
      if (!labels.empty()) {
        hasLabels = true;
        // If first label is Frame assume it is the index column
        if (labels[0] == "Frame" && indexcol_ == -1)
          indexcol_ = 0;
      }
      linebuffer = buffer.Line();
      ptr = linebuffer;
      if (ptr == 0) {
        mprinterr("Error: No data found in file.\n");
        return 1;
      }
    } else 
      // Not a recognized comment character, assume data.
      isCommentLine = false;
  }
  // Column user args start from 1
  if (indexcol_ > -1)
    mprintf("\tUsing column %i as index column.\n", indexcol_ + 1);

  // Should be at first data line. Tokenize the line.
  ntoken = buffer.TokenizeLine( SEPARATORS );
  // If # of data columns does not match # labels, clear labels.
  if ( !labels.empty() && ntoken != labels.Nargs() ) {
    labels.ClearList();
    hasLabels = false;
  }
  if (indexcol_ != -1 && indexcol_ >= ntoken) {
    mprinterr("Error: Specified index column %i is out of range (%i columns).\n",
              indexcol_+1, ntoken);
    return 1;
  }

  // Determine the type of data stored in each column. Assume numbers should
  // be read with double precision.
  MetaData md( dsname );
  DataSetList::DataListType inputSets;
  for (int col = 0; col != ntoken; ++col) {
    md.SetIdx( col+1 );
    if (hasLabels) md.SetLegend( labels[col] );
    std::string token( buffer.NextToken() );
    if (validInteger(token) || validDouble(token)) {
      // Number
      inputSets.push_back( new DataSet_double() );
    } else {
      // Assume string. Not allowed for index column.
      if (col == indexcol_) {
        mprintf("Warning: '%s' index column %i has string values. No indices will be read.\n", 
                  buffer.Filename().full(), indexcol_+1);
        indexcol_ = -1;
      }
      inputSets.push_back( new DataSet_string() );
    }
    inputSets.back()->SetMeta( md );
  }
  if (inputSets.empty()) {
    mprinterr("Error: No data detected.\n");
    return 1;
  }

  // Read in data
  while (linebuffer != 0) {
    if ( buffer.TokenizeLine( SEPARATORS ) != ntoken ) {
      PrintColumnError(buffer.LineNumber());
      break;
    }
    // Convert data in columns
    for (int i = 0; i < ntoken; ++i) {
      const char* token = buffer.NextToken();
      if (inputSets[i]->Type() == DataSet::DOUBLE)
        ((DataSet_double*)inputSets[i])->AddElement( atof(token) );
      else
        ((DataSet_string*)inputSets[i])->AddElement( std::string(token) );
    }
    //Ndata++;
    linebuffer = buffer.Line();
  }
  buffer.CloseFile();
   mprintf("\tDataFile %s has %i columns, %i lines.\n", buffer.Filename().full(),
           ntoken, buffer.LineNumber());

  // Create list containing only data sets.
  DataSetList::DataListType mySets;
  DataSet_double* Xptr = 0;
  for (int idx = 0; idx != (int)inputSets.size(); idx++)
    if ( idx != indexcol_ )
      mySets.push_back( inputSets[idx] );
    else
      Xptr = (DataSet_double*)inputSets[idx];
  std::string Xlabel;
  if (indexcol_ != -1 && indexcol_ < labels.Nargs())
    Xlabel = labels[indexcol_];
  if (Xptr == 0)
    datasetlist.AddOrAppendSets(Xlabel, DataSetList::Darray(), mySets);
  else {
    datasetlist.AddOrAppendSets(Xlabel, Xptr->Data(), mySets);
    delete Xptr;
  }

  return 0;
}

// DataIO_Std::Read_2D()
int DataIO_Std::Read_2D(std::string const& fname, 
                        DataSetList& datasetlist, std::string const& dsname)
{
  // Buffer file
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;
  mprintf("\tData will be read as a 2D square matrix.\n");
  // Skip comments
  const char* linebuffer = buffer.Line();
  while (linebuffer != 0 && linebuffer[0] == '#')
    linebuffer = buffer.Line();
  int ncols = -1;
  int nrows = 0;
  std::vector<double> matrixArray;
  while (linebuffer != 0) {
    int ntokens = buffer.TokenizeLine( SEPARATORS );
    if (ncols < 0) {
      ncols = ntokens;
      if (ntokens < 1) {
        mprinterr("Error: Could not tokenize line.\n");
        return 1;
      }
    } else if (ncols != ntokens) {
      mprinterr("Error: In 2D file, number of columns changes from %i to %i at line %i\n",
                ncols, ntokens, buffer.LineNumber());
      return 1;
    }
    for (int i = 0; i < ntokens; i++)
      matrixArray.push_back( atof( buffer.NextToken() ) );
    nrows++;
    linebuffer = buffer.Line();
  }
  if (ncols < 0) {
    mprinterr("Error: No data detected in %s\n", buffer.Filename().full());
    return 1;
  }
  if ( DetermineMatrixType( matrixArray, nrows, ncols, datasetlist, dsname )==0 ) return 1;

  return 0;
}

// DataIO_Std::Read_3D()
int DataIO_Std::Read_3D(std::string const& fname, 
                        DataSetList& datasetlist, std::string const& dsname)
{
  return 1;
}

// DataIO_Std::Read_Vector()
int DataIO_Std::Read_Vector(std::string const& fname, 
                            DataSetList& datasetlist, std::string const& dsname)
{
  // Buffer file
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;
  mprintf("\tAttempting to read vector data.\n");
  // Skip comments
  const char* linebuffer = buffer.Line();
  while (linebuffer != 0 && linebuffer[0] == '#')
    linebuffer = buffer.Line();
  // Determine format. Expect 3 (VXYZ), 6 (VXYZ OXYZ), or
  // 9 (VXYZ OXYZ VXYZ+OXYZ) values, optionally with indices.
  int ntokens = buffer.TokenizeLine( SEPARATORS );
  int ncols = ntokens; // Number of columns of vector data.
  int nv = 0;          // Number of columns to actually read from (3 or 6).
  bool hasIndex;
  if (ntokens < 1) {
    mprinterr("Error: Could not tokenize line.\n");
    return 1;
  }
  if (ncols == 3 || ncols == 6 || ncols == 9)
    hasIndex = false;
  else if (ncols == 4 || ncols == 7 || ncols == 10) {
    hasIndex = true;
    mprintf("Warning: Not reading vector data indices.\n");
  } else {
    mprinterr("Error: Expected 3, 6, or 9 columns of vector data, got %i.\n", ncols);
    return 1;
  }
  if (ncols >= 6) {
    nv = 6;
    mprintf("\tReading vector X Y Z and origin X Y Z values.\n");
  } else {
    nv = 3;
    mprintf("\tReading vector X Y Z values.\n");
  }
  // Create set
  DataSet_Vector* ds = new DataSet_Vector();
  if (ds == 0) return 1;
  ds->SetMeta( dsname );
  // Read vector data
  double vec[6];
  std::fill(vec, vec+6, 0.0);
  size_t ndata = 0;
  while (linebuffer != 0) {
    if (hasIndex)
      ntokens = sscanf(linebuffer, "%*f %lf %lf %lf %lf %lf %lf",
                       vec, vec+1, vec+2, vec+3, vec+4, vec+5);
    else
      ntokens = sscanf(linebuffer, "%lf %lf %lf %lf %lf %lf",
                       vec, vec+1, vec+2, vec+3, vec+4, vec+5);
    if (ntokens != nv) {
      mprinterr("Error: In vector file, line %i: expected %i values, got %i\n",
                buffer.LineNumber(), nv, ntokens);
      break;
    }
    ds->Add( ndata++, vec ); 
    linebuffer = buffer.Line();
  }
  return (datasetlist.AddOrAppendSets("", DataSetList::Darray(), DataSetList::DataListType(1, ds)));
}

// DataIO_Std::Read_Mat3x3()
int DataIO_Std::Read_Mat3x3(std::string const& fname, 
                            DataSetList& datasetlist, std::string const& dsname)
{
  // Buffer file
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;
  mprintf("\tAttempting to read 3x3 matrix data.\n");
  // Skip comments
  const char* linebuffer = buffer.Line();
  while (linebuffer != 0 && linebuffer[0] == '#')
    linebuffer = buffer.Line();
  // Check that number of columns (9) is correct.
  int ntokens = buffer.TokenizeLine( SEPARATORS );
  if (ntokens < 1) {
    mprinterr("Error: Could not tokenize line.\n");
    return 1;
  }
  int ncols = ntokens;
  bool hasIndex;
  if (ncols == 9)
    hasIndex = false;
  else if (ncols == 10) {
    hasIndex = true;
    mprintf("Warning: Not reading 3x3 matrix data indices.\n");
  } else {
    mprinterr("Error: Expected 9 columns of 3x3 matrix data, got %i.\n", ncols);
    return 1;
  }
  // Create data set
  DataSet_Mat3x3* ds = new DataSet_Mat3x3();
  if (ds == 0) return 1;
  ds->SetMeta( dsname );
  // Read 3x3 matrix data
  double mat[9];
  std::fill(mat, mat, 0.0);
  size_t ndata = 0;
  while (linebuffer != 0) {
    if (hasIndex)
      ntokens = sscanf(linebuffer, "%*f %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                       mat, mat+1, mat+2, mat+3, mat+4, mat+5, mat+6, mat+7, mat+8);
    else
      ntokens = sscanf(linebuffer, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                       mat, mat+1, mat+2, mat+3, mat+4, mat+5, mat+6, mat+7, mat+8);
    if (ntokens != 9) {
      mprinterr("Error: In 3x3 matrix file, line %i: expected 9 values, got %i\n",
                buffer.LineNumber(), ntokens);
      break;
    }
    ds->Add( ndata++, mat ); 
    linebuffer = buffer.Line();
  }
  return (datasetlist.AddOrAppendSets("", DataSetList::Darray(), DataSetList::DataListType(1, ds)));
}

// -----------------------------------------------------------------------------
void DataIO_Std::WriteHelp() {
  mprintf("\tinvert:     Flip X/Y axes.\n"
          "\tnoxcol:     Do not print X (index) column.\n"
          "\tnoheader:   Do not print header line.\n"
          "\tsquare2d:   Write 2D data sets in matrix-like format.\n"
          "\tnosquare2d: Write 2D data sets as '<X> <Y> <Value>'.\n");
}

// DataIO_Std::processWriteArgs()
int DataIO_Std::processWriteArgs(ArgList &argIn) {
  if (!isInverted_ && argIn.hasKey("invert"))
    isInverted_ = true;
  if (hasXcolumn_ && argIn.hasKey("noxcol"))
    hasXcolumn_ = false;
  if (writeHeader_ && argIn.hasKey("noheader"))
    writeHeader_ = false;
  bool original_square2d = square2d_;
  if (argIn.hasKey("square2d"))
    square2d_ = true;
  else if (argIn.hasKey("nosquare2d"))
    square2d_ = false;
  else
    square2d_ = original_square2d;
  return 0;
}

// WriteNameToBuffer()
void DataIO_Std::WriteNameToBuffer(CpptrajFile& fileIn, std::string const& label,
                                   int width,  bool isLeftCol)
{
  std::string temp_name = label;
  // If left aligning, add '#' to name; 
  if (isLeftCol) {
    if (temp_name[0]!='#') {
      temp_name.insert(0,"#");
      // Ensure that name will not be larger than column width.
      if ((int)temp_name.size() > width)
        temp_name.resize( width );
    }
  }
  // Replace any spaces with underscores
  for (std::string::iterator tc = temp_name.begin(); tc != temp_name.end(); ++tc)
    if ( *tc == ' ' )
      *tc = '_';
  if (width >= (int)CpptrajFile::BUF_SIZE)
    // Protect against CpptrajFile buffer overflow
    fileIn.Write(temp_name.c_str(), temp_name.size());
  else {
    // Set up header format string
    TextFormat::AlignType align;
    if (isLeftCol)
      align = TextFormat::LEFT;
    else
      align = TextFormat::RIGHT;
    TextFormat header_format(TextFormat::STRING, width, align);
    fileIn.Printf(header_format.fmt(), temp_name.c_str()); 
  }
}

// DataIO_Std::WriteData()
int DataIO_Std::WriteData(FileName const& fname, DataSetList const& SetList)
{
  int err = 0;
  if (!SetList.empty()) {
    // Open output file.
    CpptrajFile file;
    if (file.OpenWrite( fname )) return 1;
    // Base write type off first data set dimension FIXME
    if (SetList[0]->Group() == DataSet::CLUSTERMATRIX) {
      // Special case of 2D - may have sieved frames.
      err = WriteCmatrix(file, SetList);
    } else if (SetList[0]->Ndim() == 1) {
      if (isInverted_)
        err = WriteDataInverted(file, SetList);
      else
        err = WriteDataNormal(file, SetList);
    } else if (SetList[0]->Ndim() == 2)
      err = WriteData2D(file, SetList);
    else if (SetList[0]->Ndim() == 3)
      err = WriteData3D(file, SetList);
    file.CloseFile();
  }
  return err;
}

// DataIO_Std::WriteCmatrix()
int DataIO_Std::WriteCmatrix(CpptrajFile& file, DataSetList const& Sets) {
  for (DataSetList::const_iterator ds = Sets.begin(); ds != Sets.end(); ++ds)
  {
    if ( (*ds)->Group() != DataSet::CLUSTERMATRIX) {
      mprinterr("Error: Write of cluster matrix and other sets to same file not supported.\n"
                "Error: Skipping '%s'\n", (*ds)->legend());
      continue;
    }
    DataSet_Cmatrix const& cm = static_cast<DataSet_Cmatrix const&>( *(*ds) );
    int nrows = cm.Nrows();
    int col_width = std::min(3, DigitWidth( nrows ) + 1);
    int dat_width = std::max(cm.Format().Width(), (int)cm.Meta().Legend().size()) + 1;
    WriteNameToBuffer(file, "F1",               col_width, true);
    WriteNameToBuffer(file, "F2",               col_width, false);
    WriteNameToBuffer(file, cm.Meta().Legend(), dat_width, false);
    file.Printf("\n");
    TextFormat col_fmt(TextFormat::INTEGER, col_width);
    TextFormat dat_fmt = cm.Format();
    dat_fmt.SetFormatAlign(TextFormat::RIGHT);
    dat_fmt.SetFormatWidth( dat_width );
    std::string total_fmt = col_fmt.Fmt() + col_fmt.Fmt() + dat_fmt.Fmt() + "\n";
    //mprintf("DEBUG: format '%s'\n", total_fmt.c_str());
    ClusterSieve::SievedFrames const& frames = cm.FramesToCluster();
    int ntotal = (int)frames.size();
    for (int idx1 = 0; idx1 != ntotal; idx1++) {
      int row = frames[idx1];
      for (int idx2 = idx1 + 1; idx2 != ntotal; idx2++) {
        int col = frames[idx2];
        file.Printf(total_fmt.c_str(), row+1, col+1, cm.GetFdist(col, row)); 
      }
    }
  }
  return 0;
}

// DataIO_Std::WriteDataNormal()
int DataIO_Std::WriteDataNormal(CpptrajFile& file, DataSetList const& Sets) {
  // Assume all 1D data sets.
  if (Sets.empty() || CheckAllDims(Sets, 1)) return 1;
  // For this output to work the X-dimension of all sets needs to match.
  // The most important things for output are min and step so just check that.
  // Use X dimension of set 0 for all set dimensions.
  CheckXDimension( Sets );
  // TODO: Check for empty dim.
  DataSet* Xdata = Sets[0];
  Dimension const& Xdim = static_cast<Dimension const&>( Xdata->Dim(0) );
  int xcol_width = Xdim.Label().size() + 1; // Only used if hasXcolumn_
  if (xcol_width < 8) xcol_width = 8;
  int xcol_precision = 3;

  // Determine size of largest DataSet.
  size_t maxFrames = DetermineMax( Sets );

  // Set up X column.
  TextFormat x_col_format(XcolFmt());
  if (hasXcolumn_) {
    if (XcolPrecSet()) {
      x_col_format = TextFormat(XcolFmt(), XcolWidth(), XcolPrec());
    } else {
      // Create format string for X column based on dimension in first data set.
      // Adjust X col precision as follows: if the step is set and has a 
      // fractional component set the X col width/precision to either the data
      // width/precision or the current width/precision, whichever is larger. If
      // the set is XYMESH but step has not been set (so we do not know spacing 
      // between X values) use default precision. Otherwise the step has no
      // fractional component so make the precision zero.
      double step_i;
      double step_f = modf( Xdim.Step(), &step_i );
      double min_f  = modf( Xdim.Min(),  &step_i );
      if (Xdim.Step() > 0.0 && (step_f > 0.0 || min_f > 0.0)) {
        xcol_precision = std::max(xcol_precision, Xdata->Format().Precision());
        xcol_width = std::max(xcol_width, Xdata->Format().Width());
      } else if (Xdata->Type() != DataSet::XYMESH)
        xcol_precision = 0;
      x_col_format.SetCoordFormat( maxFrames, Xdim.Min(), Xdim.Step(), xcol_width, xcol_precision );
    }
  } else {
    // If not writing an X-column, no leading space for the first dataset.
    Xdata->SetupFormat().SetFormatAlign( TextFormat::RIGHT );
  }

  // Write header to buffer
  std::vector<int> Original_Column_Widths;
  if (writeHeader_) {
    // If x-column present, write x-label
    if (hasXcolumn_)
      WriteNameToBuffer( file, Xdim.Label(), xcol_width, true );
    // To prevent truncation of DataSet legends, adjust the width of each
    // DataSet if necessary.
    for (DataSetList::const_iterator ds = Sets.begin(); ds != Sets.end(); ++ds) {
      // Record original column widths in case they are changed.
      Original_Column_Widths.push_back( (*ds)->Format().Width() );
      int colLabelSize;
      if (ds == Sets.begin() && !hasXcolumn_)
        colLabelSize = (int)(*ds)->Meta().Legend().size() + 1;
      else
        colLabelSize = (int)(*ds)->Meta().Legend().size();
      //mprintf("DEBUG: Set '%s', fmt width= %i, colWidth= %i, colLabelSize= %i\n",
      //        (*ds)->legend(), (*ds)->Format().Width(), (*ds)->Format().ColumnWidth(),
      //        colLabelSize);
      if (colLabelSize >= (*ds)->Format().ColumnWidth())
        (*ds)->SetupFormat().SetFormatWidth( colLabelSize );
    }
    // Write dataset names to header, left-aligning first set if no X-column
    DataSetList::const_iterator set = Sets.begin();
    if (!hasXcolumn_)
      WriteNameToBuffer( file, (*set)->Meta().Legend(), (*set)->Format().ColumnWidth(), true  );
    else
      WriteNameToBuffer( file, (*set)->Meta().Legend(), (*set)->Format().ColumnWidth(), false );
    ++set;
    for (; set != Sets.end(); ++set) 
      WriteNameToBuffer( file, (*set)->Meta().Legend(), (*set)->Format().ColumnWidth(), false );
    file.Printf("\n"); 
  }

  // Write Data
  DataSet::SizeArray positions(1);
  for (positions[0] = 0; positions[0] < maxFrames; positions[0]++) {
    // Output Frame for each set
    if (hasXcolumn_)
      file.Printf( x_col_format.fmt(), Xdata->Coord(0, positions[0]) );
    for (DataSetList::const_iterator set = Sets.begin(); set != Sets.end(); ++set) 
      (*set)->WriteBuffer(file, positions);
    file.Printf("\n"); 
  }
  // Restore original column widths if necessary
  if (!Original_Column_Widths.empty())
    for (unsigned int i = 0; i != Original_Column_Widths.size(); i++)
      Sets[i]->SetupFormat().SetFormatWidth( Original_Column_Widths[i] );
  return 0;
}

// DataIO_Std::WriteDataInverted()
int DataIO_Std::WriteDataInverted(CpptrajFile& file, DataSetList const& Sets)
{
  if (Sets.empty() || CheckAllDims(Sets, 1)) return 1;
  // Determine size of largest DataSet.
  size_t maxFrames = DetermineMax( Sets );
  // Write each set to a line.
  DataSet::SizeArray positions(1);
  for (DataSetList::const_iterator set = Sets.begin(); set != Sets.end(); ++set) {
    // Write dataset name as first column.
    WriteNameToBuffer( file, (*set)->Meta().Legend(), (*set)->Format().ColumnWidth(), false); 
    // Write each frame to subsequent columns
    for (positions[0] = 0; positions[0] < maxFrames; positions[0]++) 
      (*set)->WriteBuffer(file, positions);
    file.Printf("\n");
  }
  return 0;
}

// DataIO_Std::WriteData2D()
int DataIO_Std::WriteData2D( CpptrajFile& file, DataSetList const& setList) 
{
  int err = 0;
  for (DataSetList::const_iterator set = setList.begin(); set != setList.end(); ++set)
  {
    if (set != setList.begin()) file.Printf("\n");
    err += WriteSet2D( *(*set), file );
  }
  return err;
}

// DataIO_Std::WriteSet2D()
int DataIO_Std::WriteSet2D( DataSet const& setIn, CpptrajFile& file ) {
  if (setIn.Ndim() != 2) {
    mprinterr("Internal Error: DataSet %s in DataFile %s has %zu dimensions, expected 2.\n",
              setIn.legend(), file.Filename().full(), setIn.Ndim());
    return 1;
  }
  DataSet_2D const& set = static_cast<DataSet_2D const&>( setIn );
  int xcol_width = 8;
  int xcol_precision = 3;
  Dimension const& Xdim = static_cast<Dimension const&>(set.Dim(0));
  Dimension const& Ydim = static_cast<Dimension const&>(set.Dim(1));
  if (Xdim.Step() == 1.0) xcol_precision = 0;
  
  DataSet::SizeArray positions(2);
  TextFormat ycoord_fmt(XcolFmt()), xcoord_fmt(XcolFmt());
  if (square2d_) {
    // Print XY values in a grid:
    // x0y0 x1y0 x2y0
    // x0y1 x1y1 x2y1
    // x0y2 x1y2 x2y2
    // If file has header, top-left value will be '#<Xlabel>-<Ylabel>',
    // followed by X coordinate values.
    if (writeHeader_) {
      ycoord_fmt.SetCoordFormat( set.Nrows(), Ydim.Min(), Ydim.Step(), xcol_width, xcol_precision );
      std::string header;
      if (Xdim.Label().empty() && Ydim.Label().empty())
        header = "#Frame";
      else
        header = "#" + Xdim.Label() + "-" + Ydim.Label();
      WriteNameToBuffer( file, header, xcol_width, true );
      xcoord_fmt.SetCoordFormat( set.Ncols(), Xdim.Min(), Xdim.Step(),
                                 set.Format().ColumnWidth(), xcol_precision );
      for (size_t ix = 0; ix < set.Ncols(); ix++)
        file.Printf( xcoord_fmt.fmt(), set.Coord(0, ix) );
      file.Printf("\n");
    }
    for (positions[1] = 0; positions[1] < set.Nrows(); positions[1]++) {
      if (writeHeader_)
        file.Printf( ycoord_fmt.fmt(), set.Coord(1, positions[1]) );
      for (positions[0] = 0; positions[0] < set.Ncols(); positions[0]++)
        set.WriteBuffer( file, positions );
      file.Printf("\n");
    }
  } else {
    // Print X Y Values
    // x y val(x,y)
    if (writeHeader_)
      file.Printf("#%s %s %s\n", Xdim.Label().c_str(), 
                  Ydim.Label().c_str(), set.legend());
    if (XcolPrecSet()) {
      xcoord_fmt = TextFormat(XcolFmt(), XcolWidth(), XcolPrec());
      ycoord_fmt = xcoord_fmt;
    } else {
      xcoord_fmt.SetCoordFormat( set.Ncols(), Xdim.Min(), Xdim.Step(), 8, 3 );
      ycoord_fmt.SetCoordFormat( set.Nrows(), Ydim.Min(), Ydim.Step(), 8, 3 );
    }
    std::string xy_fmt = xcoord_fmt.Fmt() + " " + ycoord_fmt.Fmt() + " ";
    for (positions[1] = 0; positions[1] < set.Nrows(); ++positions[1]) {
      for (positions[0] = 0; positions[0] < set.Ncols(); ++positions[0]) {
        file.Printf( xy_fmt.c_str(), set.Coord(0, positions[0]), set.Coord(1, positions[1]) );
        set.WriteBuffer( file, positions );
        file.Printf("\n");
      }
    }
  }
  return 0;
}

// DataIO_Std::WriteData3D()
int DataIO_Std::WriteData3D( CpptrajFile& file, DataSetList const& setList) 
{
  int err = 0;
  for (DataSetList::const_iterator set = setList.begin(); set != setList.end(); ++set)
  {
    if (set != setList.begin()) file.Printf("\n");
    err += WriteSet3D( *(*set), file );
  }
  return err;
}

// DataIO_Std::WriteSet3D()
int DataIO_Std::WriteSet3D( DataSet const& setIn, CpptrajFile& file ) {
  if (setIn.Ndim() != 3) {
    mprinterr("Internal Error: DataSet %s in DataFile %s has %zu dimensions, expected 3.\n",
              setIn.legend(), file.Filename().full(), setIn.Ndim());
    return 1;
  }
  DataSet_3D const& set = static_cast<DataSet_3D const&>( setIn );
  Dimension const& Xdim = static_cast<Dimension const&>(set.Dim(0));
  Dimension const& Ydim = static_cast<Dimension const&>(set.Dim(1));
  Dimension const& Zdim = static_cast<Dimension const&>(set.Dim(2));
  //if (Xdim.Step() == 1.0) xcol_precision = 0;
  
  // Print X Y Z Values
  // x y z val(x,y,z)
  DataSet::SizeArray pos(3);
  if (writeHeader_)
    file.Printf("#%s %s %s %s\n", Xdim.Label().c_str(), 
                Ydim.Label().c_str(), Zdim.Label().c_str(), set.legend());
  std::string xyz_fmt;
  if (XcolPrecSet()) {
    TextFormat nfmt( XcolFmt(), XcolWidth(), XcolPrec() );
    xyz_fmt = nfmt.Fmt() + " " + nfmt.Fmt() + " " + nfmt.Fmt() + " ";
  } else {
    TextFormat xfmt( XcolFmt(), set.NX(), Xdim.Min(), Xdim.Step(), 8, 3 );
    TextFormat yfmt( XcolFmt(), set.NY(), Ydim.Min(), Ydim.Step(), 8, 3 );
    TextFormat zfmt( XcolFmt(), set.NZ(), Zdim.Min(), Zdim.Step(), 8, 3 );
    xyz_fmt = xfmt.Fmt() + " " + yfmt.Fmt() + " " + zfmt.Fmt() + " ";
  }
  for (pos[2] = 0; pos[2] < set.NZ(); ++pos[2]) {
    for (pos[1] = 0; pos[1] < set.NY(); ++pos[1]) {
      for (pos[0] = 0; pos[0] < set.NX(); ++pos[0]) {
        file.Printf( xyz_fmt.c_str(), set.Coord(0, pos[0]),
                     set.Coord(1, pos[1]), set.Coord(2, pos[2]) );
        set.WriteBuffer( file, pos );
        file.Printf("\n");
      }
    }
  }
  return 0;
}

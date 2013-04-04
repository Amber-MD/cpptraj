#include <cstdlib> // atoi, atof
#include <cstring> // strchr
#include <cctype>  // isdigit, isalpha
#include "DataIO_Std.h"
#include "CpptrajStdio.h" 
#include "StringRoutines.h" // SetStringFormatString
#include "BufferedLine.h"

// CONSTRUCTOR
DataIO_Std::DataIO_Std() : 
  writeHeader_(true),
  square2d_(false),
  ymin_(1),
  ystep_(1)
{}

static void PrintColumnError() {
  mprinterr("Error: Number of columns in file changes.\n");
}

// DataIO_Std::ReadData()
int DataIO_Std::ReadData(std::string const& fname, DataSetList& datasetlist) {
  ArgList labels;
  bool hasLabels = false;
  std::vector<DataSet*> DsetList;
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
    DataSet* dset = 0;
    if ( isalpha( token[0] ) ) 
    {
      //mprintf(" STRING!\n");
      // STRING columns cannot be index columns
      if ( col == indexcol ) {
        mprinterr("Error: DataFile %s index column %i has string values.\n", 
                  buffer.Filename().full(), indexcol+1);
        return 1;
      }
      dset = datasetlist.AddSetIdx( DataSet::STRING, buffer.Filename().Base(), col+1 );
    } else if ( isdigit( token[0] ) || 
                token[0]=='+' || 
                token[0]=='-' ||
                token[0]=='.'   )
    {
      if ( strchr( token, '.' ) != 0 ) {
        //mprintf(" DOUBLE!\n");
        if ( col != indexcol )
          dset = datasetlist.AddSetIdx( DataSet::DOUBLE, buffer.Filename().Base(), col+1 );
        //else
        //  indextype = DataSet::DOUBLE;
      } else {
        //mprintf(" INTEGER!\n");
        if (col != indexcol)
          dset = datasetlist.AddSetIdx( DataSet::INT, buffer.Filename().Base(), col+1 );
        //else
        //  indextype = DataSet::INT;
      }
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
        case DataSet::INT: 
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
  std::string ylabel = argIn.GetStringKey("ylabel");
  if (!ylabel.empty()) y_label_.assign(ylabel);
  ymin_ = argIn.getKeyDouble("ymin",ymin_);
  ystep_ = argIn.getKeyDouble("ystep",ystep_);
 
  if (argIn.hasKey("noheader"))
    writeHeader_ = false;
  if (argIn.hasKey("square2d"))
    square2d_ = true;

  return 0;
}

// DataIO_Std::WriteNameToBuffer()
void DataIO_Std::WriteNameToBuffer(CpptrajFile& fileIn, DataSet* DS, bool leftAlign) 
{
  std::string temp_name = DS->Legend();
  // If left aligning, add '#' to name; ensure that name will not be
  // larger than column width.
  if (leftAlign) {
    if (temp_name[0]!='#')
      temp_name.insert(0,"#");
  }
  int width = DS->Width();
  if ((int)temp_name.size() > width)
    temp_name.resize( width );
  // Replace any spaces with underscores
  for (std::string::iterator tc = temp_name.begin(); tc != temp_name.end(); ++tc)
    if ( *tc == ' ' )
      *tc = '_';
  // Set up header format string
  std::string header_format;
  SetStringFormatString(header_format, width, leftAlign);
  fileIn.Printf(header_format.c_str(), temp_name.c_str());
}

// DataIO_Std::WriteData()
int DataIO_Std::WriteData(std::string const& fname, DataSetList &SetList) {
  std::string x_header_fmt;
  DataSetList::const_iterator set;

  CpptrajFile file;
  if (file.OpenWrite( fname )) return 1;
  // Create format string for X column. If xstep is 1 set precision to 0
  // NOTE: only is hasXcolumn?
  if (xstep_ == 1) xcol_precision_ = 0;
  SetupXcolumn();

  // If not writing an X-column, set the format for the first dataset
  // to left-aligned.
  if (!hasXcolumn_) {
    set = SetList.begin();
    (*set)->SetDataSetFormat( true );
  }

  // Write header to buffer
  if (writeHeader_) {
    // If x-column present, write x-label
    if (hasXcolumn_) {
      // Insert leading '#'character. Ensure the result is no greater 
      // than xcol_width.
      x_label_.insert(0,"#");
      x_label_.resize( xcol_width_, ' ');
      SetStringFormatString(x_header_fmt, xcol_width_, true);
      file.Printf(x_header_fmt.c_str(), x_label_.c_str());
    }
    // Write dataset names to header, left-aligning first set if no X-column
    set = SetList.begin();
    if (!hasXcolumn_)
      WriteNameToBuffer( file, *set, true );
    else
      WriteNameToBuffer( file, *set, false );
    ++set;
    for (; set != SetList.end(); ++set) 
      WriteNameToBuffer( file, *set, false );
    file.Printf("\n"); 
  }

  // Write Data
  for (int frame=0; frame < maxFrames_; frame++) {
    // If not printing empty frames, make sure that every set has data
    // at this frame.
    if (!printEmptyFrames_) {
      bool emptyFrames = false;
      for (set = SetList.begin(); set != SetList.end(); set++) {
        if ( (*set)->FrameIsEmpty(frame) ) {
          emptyFrames = true;
          break;
        }
      }
      if (emptyFrames) continue;
    }
    // Output Frame
    if (hasXcolumn_)
      PrintX(file, frame); 
    for (set = SetList.begin(); set != SetList.end(); set++) 
      (*set)->WriteBuffer(file,frame);
    file.Printf("\n"); 
  }
  file.CloseFile();
  return 0;
}

// DataIO_Std::WriteDataInverted()
int DataIO_Std::WriteDataInverted(std::string const& fname, DataSetList &SetList) {
  DataSetList::const_iterator set;
  std::string dset_name, x_header_fmt;
  CpptrajFile file;
  if (file.OpenWrite( fname )) return 1;

  for (set = SetList.begin(); set != SetList.end(); set++) {
    // if specified check for empty frames in the set
    if (!printEmptyFrames_) {
      bool emptyFrames = false; 
      for (int frame=0; frame < maxFrames_; frame++) {
        if ((*set)->FrameIsEmpty(frame) ) {
          emptyFrames = true; 
          break;
        }
      }
      if (emptyFrames) continue;
    }
    // Write dataset name as first column.
    WriteNameToBuffer( file, *set, false); 
    // Write each frame to subsequent columns
    for (int frame=0; frame<maxFrames_; frame++) 
      (*set)->WriteBuffer(file,frame);
    file.Printf("\n");
  }
  file.CloseFile();
  return 0;
}

int DataIO_Std::WriteData2D( std::string const& fname, DataSet& set ) {
  std::vector<int> dimensions; 
  CpptrajFile file;
  if (file.OpenWrite( fname )) return 1;
  set.GetDimensions(dimensions);
  if (dimensions.size() != 2) {
    mprinterr("Internal Error: DataSet %s in DataFile %s has %zu dimensions, expected 2.\n",
              set.Legend().c_str(), file.Filename().full(), dimensions.size());
    return 1;
  }
  
  if (square2d_) {
    // Print XY values in a grid
    // x0y0 x0y1 x0y2
    // x1y0 x1y1 x1y2
    std::string headerstring;
    if (writeHeader_) {
      SetIntegerFormatString( headerstring, set.Width(), false);
      if (hasXcolumn_)
        file.Printf("%-*s ",set.Width(), "#Frame");
      for (int iy = 0; iy < dimensions[1]; ++iy)
        file.Printf(headerstring.c_str(), iy+OUTPUTFRAMESHIFT);
      file.Printf("\n");
    }
    for (int ix = 0; ix < dimensions[0]; ++ix) {
      if (hasXcolumn_)
        file.Printf(headerstring.c_str(), ix+OUTPUTFRAMESHIFT);
      for (int iy = 0; iy < dimensions[1]; ++iy) {
        set.Write2D( file, ix, iy);
      }
      file.Printf("\n");
    }
  } else {
    // Print X Y Value
    // Print dataset name
    file.Printf("#%s\n", set.Legend().c_str());
  
    for (int ix = 0; ix < dimensions[0]; ++ix) {
      double xcoord = (xstep_ * (double)ix) + xmin_;
      for (int iy = 0; iy < dimensions[1]; ++iy) {
        double ycoord = (ystep_ * (double)iy) + ymin_;
        file.Printf("%8.3f %8.3f", xcoord, ycoord);
        set.Write2D( file, ix, iy );
        file.Printf("\n");
      }
    }
  }
  return 0;
}

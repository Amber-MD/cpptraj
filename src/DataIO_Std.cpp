#include "DataIO_Std.h"
#include "CpptrajStdio.h" // SetStringFormatString
// Eventually fold this into CpptrajFile?
#include "CharBuffer.h"

// CONSTRUCTOR
DataIO_Std::DataIO_Std() : 
  writeHeader_(true),
  square2d_(false),
  ymin_(1),
  ystep_(1)
{}

// DataIO_Std::processWriteArgs()
int DataIO_Std::processWriteArgs(ArgList &argIn) {
  char *ylabel = argIn.getKeyString("ylabel",NULL);
  if (ylabel!=NULL) y_label_.assign(ylabel);
  ymin_ = argIn.getKeyDouble("ymin",ymin_);
  ystep_ = argIn.getKeyDouble("ystep",ystep_);
 
  if (argIn.hasKey("noheader"))
    writeHeader_ = false;
  if (argIn.hasKey("square2d"))
    square2d_ = true;

  return 0;
}

// DataIO_Std::WriteNameToBuffer()
void DataIO_Std::WriteNameToBuffer(CharBuffer& cbuffer, DataSet* DS, bool leftAlign) 
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
  // Set up header format string
  std::string header_format;
  SetStringFormatString(header_format, width, leftAlign);
  cbuffer.Sprintf(header_format.c_str(), temp_name.c_str());
}

// DataIO_Std::WriteData()
int DataIO_Std::WriteData(DataSetList &SetList) {
  CharBuffer buffer;
  std::string x_header_fmt;
  DataSetList::const_iterator set;

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
      buffer.Sprintf(x_header_fmt.c_str(), x_label_.c_str());
    }
    // Write dataset names to header, left-aligning first set if no X-column
    set = SetList.begin();
    if (!hasXcolumn_)
      WriteNameToBuffer( buffer, *set, true );
    else
      WriteNameToBuffer( buffer, *set, false );
    ++set;
    for (; set != SetList.end(); ++set) 
      WriteNameToBuffer( buffer, *set, false );
    buffer.NewLine();
  }

  // Ensure buffer has enough space for all data in dataset list
  // Each dataset takes up (maxFrames * (width + space + newline))
  size_t total_datasize = 0;
  for (set = SetList.begin(); set != SetList.end(); set++) 
    total_datasize += (size_t)( maxFrames_ * ((*set)->Width() + 2) );
  // If X-column present, xcol_width per frame
  if (hasXcolumn_)
    total_datasize += (maxFrames_ * xcol_width_);
  buffer.Reallocate( total_datasize );

  // Write Data to buffer
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
    if (hasXcolumn_) {
      double xcoord = (xstep_ * (double)frame) + xmin_;
      buffer.WriteDouble(x_format_.c_str(), xcoord);
    }
    for (set = SetList.begin(); set != SetList.end(); set++) 
      (*set)->WriteBuffer(buffer,frame);
    buffer.NewLine();
  }
  //buffer.DumpBuffer();
  // Write buffer to file
  IO->Write(buffer.Buffer(),buffer.CurrentSize());
  return 0;
}

// DataIO_Std::WriteDataInverted()
int DataIO_Std::WriteDataInverted(DataSetList &SetList) {
  CharBuffer buffer;
  DataSetList::const_iterator set;
  std::string dset_name, x_header_fmt;

  // Determine max output file size
  // Each set is (number of frames * width) + newline))
  size_t dataFileSize = 0;
  for (set = SetList.begin(); set != SetList.end(); set++) {
    dataFileSize += (size_t)((maxFrames_ * (*set)->Width()) + 1);
    // Add the size of the dataset name, will be x column
    dataFileSize += (size_t)(*set)->Width();
    //dset_name.assign( (*set)->Name() );
    //size_t dset_name_size = dset_name.size();
    //if ((size_t)(*set)->Width() > dset_name_size)
    //  dset_name_size = (*set)->Width();
    //dataFileSize += dset_name_size;
  } 
  buffer.Allocate( dataFileSize );

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
    WriteNameToBuffer(buffer, *set, false); 
    // Write each frame to subsequent columns
    for (int frame=0; frame<maxFrames_; frame++) 
      (*set)->WriteBuffer(buffer,frame);
    buffer.NewLine();
  }
  IO->Write(buffer.Buffer(),buffer.CurrentSize());
  return 0;
}

int DataIO_Std::WriteData2D( DataSet& set ) {
  std::vector<int> dimensions; 
  const char* headerfmt; 

  set.GetDimensions(dimensions);
  if (dimensions.size() != 2) {
    mprinterr("Internal Error: DataSet %s in DataFile %s has %zu dimensions, expected 2.\n",
              set.c_str(), Name(), dimensions.size());
    return 1;
  }
  

  if (square2d_) {
    // Print XY values in a grid
    // x0y0 x0y1 x0y2
    // x1y0 x1y1 x1y2
    if (writeHeader_) {
      std::string headerstring;
      SetIntegerFormatString( headerstring, set.Width(), false);
      headerfmt = headerstring.c_str();
      Printf("%-*s ",set.Width(), "#Frame");
      for (int iy = 0; iy < dimensions[1]; ++iy)
        Printf(headerfmt, iy+OUTPUTFRAMESHIFT);
      Printf("\n");
    }
    for (int ix = 0; ix < dimensions[0]; ++ix) {
      if (writeHeader_)
        Printf(headerfmt, ix+OUTPUTFRAMESHIFT);
      for (int iy = 0; iy < dimensions[1]; ++iy) {
        set.Write2D( *this, ix, iy);
      }
      Printf("\n");
    }
  } else {
    // Print X Y Value
    // Print dataset name
    Printf("#%s\n", set.c_str());
  
    for (int ix = 0; ix < dimensions[0]; ++ix) {
      double xcoord = (xstep_ * (double)ix) + xmin_;
      for (int iy = 0; iy < dimensions[1]; ++iy) {
        double ycoord = (ystep_ * (double)iy) + ymin_;
        Printf("%8.3f %8.3f", xcoord, ycoord);
        set.Write2D( *this, ix, iy );
        Printf("\n");
      }
    }
  }
  return 0;
}


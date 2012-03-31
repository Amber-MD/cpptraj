#include "DataIO_Std.h"
#include "CpptrajStdio.h" // SetStringFormatString
// Eventually fold this into CpptrajFile?
#include "CharBuffer.h"

// CONSTRUCTOR
StdDataFile::StdDataFile() {
  writeHeader_ = true;
}

// StdDataFile::processWriteArgs()
int StdDataFile::processWriteArgs(ArgList &argIn) {
  if (argIn.hasKey("noheader"))
    writeHeader_ = false;
  return 0;
}

// StdDataFile::WriteData()
int StdDataFile::WriteData(DataSetList &SetList) {
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
    // Write dataset names to header
    for (set = SetList.begin(); set != SetList.end(); set++) 
      (*set)->WriteNameToBuffer( buffer );
    buffer.NewLine();
  }

  // Ensure buffer has enough space for all data in dataset list
  // Each dataset takes up (maxFrames * (width + newline))
  size_t total_datasize = 0;
  for (set = SetList.begin(); set != SetList.end(); set++) 
    total_datasize += (size_t)( maxFrames_ * ((*set)->Width() + 1) );
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

// StdDataFile::WriteDataInverted()
int StdDataFile::WriteDataInverted(DataSetList &SetList) {
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
    (*set)->WriteNameToBuffer(buffer); 
    // Write each frame to subsequent columns
    for (int frame=0; frame<maxFrames_; frame++) 
      (*set)->WriteBuffer(buffer,frame);
    buffer.NewLine();
  }
  IO->Write(buffer.Buffer(),buffer.CurrentSize());
  return 0;
}

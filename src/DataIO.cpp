#include "DataIO.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataIO::DataIO() {
  maxFrames_ = 0;

  hasXcolumn_ = true;
  xcol_width_ = 0;
  xcol_precision_ = 3;
  x_label_ = "Frame";
  xmin_ = 1;
  xstep_ = 1;

  printEmptyFrames_ = true;
}

// Copy Constructor
DataIO::DataIO(const DataIO &rhs) :
  CpptrajFile(rhs)
{
  maxFrames_ = rhs.maxFrames_;
  hasXcolumn_ = rhs.hasXcolumn_;
  xcol_width_ = rhs.xcol_width_;
  xcol_precision_ = rhs.xcol_precision_;
  x_label_ = rhs.x_label_;
  xmin_ = rhs.xmin_;
  xstep_ = rhs.xstep_;
  printEmptyFrames_ = rhs.printEmptyFrames_;
}

// Assignment
DataIO &DataIO::operator=(const DataIO &rhs) {
  // Self
  if (this == &rhs) return *this;
  // Base
  CpptrajFile::operator=(rhs);
  // Deallocate
  // Allocate and copy
  maxFrames_ = rhs.maxFrames_;
  hasXcolumn_ = rhs.hasXcolumn_;
  xcol_width_ = rhs.xcol_width_;
  xcol_precision_ = rhs.xcol_precision_;
  x_label_ = rhs.x_label_;
  xmin_ = rhs.xmin_;
  xstep_ = rhs.xstep_;
  printEmptyFrames_ = rhs.printEmptyFrames_;
  return *this;
}

// DataIO::SetDebug()
void DataIO::SetDebug(int debugIn) {
  debug_ = debugIn;
}

// DataIO::SetMaxFrames()
void DataIO::SetMaxFrames(int maxIn) {
  maxFrames_ = maxIn;
}

// DataIO::ProcessCommonArgs()
int DataIO::ProcessCommonArgs(ArgList &argIn) {
  if (argIn.hasKey("noxcol")) hasXcolumn_ = false;
  if (argIn.hasKey("noemptyframes")) printEmptyFrames_ = false;
  char *xlabel = argIn.getKeyString("xlabel",NULL);
  if (xlabel!=NULL) x_label_.assign(xlabel);
  xmin_ = argIn.getKeyDouble("xmin",xmin_);
  xstep_ = argIn.getKeyDouble("xstep",xstep_);
  return 0;
}  

// DataIO::SetupXcolumn()
/** Set the x-column format string. First check the maximum digit width
  * based on the calculated maximum x_value from maxframes, xstep, and
  * xmin. Then check how many digits the precision might add. The minimum
  * x-column width is set at 8. The x-column is left-aligned (no leading
  * space) by default.
  */
void DataIO::SetupXcolumn() {
  // Determine the character width necessary to hold the largest X value
  int max_xval = maxFrames_;
  if (xstep_ > 1)
    max_xval *= (int)xstep_;
  max_xval += (int)xmin_;
  xcol_width_ = DigitWidth( max_xval );
  // If the width for the x column plus the characters needed for precision
  // (plus 1 for decimal point) would be greater than 8, increment the 
  // X column width by (precision+1).
  if (xcol_precision_ != 0) {
    int precision_width = xcol_width_ + xcol_precision_ + 1;
    if ( precision_width > 8) xcol_width_ = precision_width;
  }
  // Default width for x col is at least 8
  if (xcol_width_ < 8) xcol_width_ = 8;
  // Set X column data format string, left-aligned (no leading space)
  SetDoubleFormatString(x_format_, xcol_width_, xcol_precision_, 0, true);
}

#include <cmath> // log10, fabs
#include "DataIO.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // DoubleFormatString

// CONSTRUCTOR
DataIO::DataIO() :
  maxFrames_(0),
  debug_(0),
  hasXcolumn_(true),
  xcol_width_(0),
  xcol_precision_(3),
  x_label_("Frame"),
  xmin_(1.0),
  xstep_(1.0),
  xoffset_(0.0),
  printEmptyFrames_(true)
{}

// Copy Constructor
DataIO::DataIO(const DataIO &rhs) : 
  maxFrames_(rhs.maxFrames_),
  debug_(rhs.debug_),
  hasXcolumn_(rhs.hasXcolumn_),
  xcol_width_(rhs.xcol_width_),
  xcol_precision_(rhs.xcol_precision_),
  x_label_(rhs.x_label_),
  xmin_(rhs.xmin_),
  xstep_(rhs.xstep_),
  xoffset_(rhs.xoffset_),
  printEmptyFrames_(rhs.printEmptyFrames_)
{}

// Assignment
DataIO &DataIO::operator=(const DataIO &rhs) {
  // Self
  if (this == &rhs) return *this;
  // Deallocate
  // Allocate and copy
  maxFrames_ = rhs.maxFrames_;
  hasXcolumn_ = rhs.hasXcolumn_;
  xcol_width_ = rhs.xcol_width_;
  xcol_precision_ = rhs.xcol_precision_;
  x_label_ = rhs.x_label_;
  xmin_ = rhs.xmin_;
  xstep_ = rhs.xstep_;
  xoffset_ = rhs.xoffset_;
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
  std::string label = argIn.GetStringKey("xlabel");
  if (!label.empty()) x_label_ = label;
  xmin_ = argIn.getKeyDouble("xmin",xmin_);
  xstep_ = argIn.getKeyDouble("xstep",xstep_);
  if (argIn.Contains("time")) {
    xstep_ = argIn.getKeyDouble("time",xstep_);
    xmin_ = 0.0;
    xoffset_ = 1.0;
  }
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
  int max_xval = maxFrames_ + (int)xoffset_;
  if (xstep_ > 1)
    max_xval *= (int)xstep_;
  max_xval += (int)xmin_;
  xcol_width_ = DigitWidth( max_xval );
  // Check if the precision is enough to support the step size
  if (xstep_ < 1.0) {
    double precision_exponent = fabs( log10( xstep_ ) );
    ++precision_exponent;
    int prec_exp_width = (int)precision_exponent; // Cast to int implicitly rounds down
    if (prec_exp_width > xcol_precision_)
      xcol_precision_ = prec_exp_width;
  }
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

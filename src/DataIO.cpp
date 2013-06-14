#include "DataIO.h"
#include "StringRoutines.h"

std::string DataIO::SetupCoordFormat(size_t maxFrames, Dimension const& dim, 
                                     int default_width, int default_precision)
{
  return SetupCoordFormat(maxFrames, dim, default_width, default_precision, true);
}

std::string DataIO::SetupCoordFormat(size_t maxFrames, Dimension const& dim, 
                                     int default_width, int default_precision,
                                     bool leftAligned) 
{
  int col_precision = default_precision;
  // Determine maximum coordinate.
  double maxCoord = (dim.Step() * (double)maxFrames) + dim.Min();
  // Determine character width necessary to hold largest coordinate.
  int col_width = DigitWidth( (long int)maxCoord );
  // Check if the precision is enough to support the step size.
  if (dim.Step() < 1.0) {
    int prec_exp_width = FloatWidth( dim.Step() );
    if (prec_exp_width > col_precision)
      col_precision = prec_exp_width;
  }
  // If the width for the column plus the characters needed for precision
  // (plus 1 for decimal point) would be greated than default_width, increment 
  // the column width by (precision+1).
  if (col_precision != 0) {
    int precision_width = col_width + col_precision + 1;
    if (precision_width > default_width) col_width = precision_width;
  }
  // Default width for column is at least default_width.
  if (col_width < default_width) col_width = default_width;
  // Set column data format string, left-aligned (no leading space).
  return SetDoubleFormatString( col_width, col_precision, 0, leftAligned );
}

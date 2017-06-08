#include "TextFormat.h"
#include "StringRoutines.h"

char TextFormat::TypeChar_[] = { 'f', 'E', 'g', 'i', 's' };

const char* TextFormat::TypeDesc_[] = {
  "DOUBLE", "SCIENTIFIC", "GENERAL", "INTEGER", "STRING"
};

// TODO benchmark - will using a big buffer and C string routines be better?
void TextFormat::SetFormatString() {
  std::string width_arg, prec_arg, left_arg, long_arg;
  // Set width and/or precision if applicable.
  if (width_ > 0)
    width_arg.assign( integerToString( width_ ) );
  if (precision_ > -1 && (int)type_ < (int)INTEGER)
    prec_arg.assign( "." + integerToString( precision_ ) );
  if (align_ == LEFT)
    left_arg.assign("-");
  if (isLong_)
    long_arg.assign("l");
  // Set format string.
  if (align_ == LEADING_SPACE) {
    fmt_.assign(" ");
    colwidth_ = 1;
  } else {
    fmt_.clear();
    colwidth_ = 0;
  }
  for (int i = 0; i != nelements_; i++) {
    if (i != 0) {
      fmt_.append(" "); // TODO: Option to not have spaces in-between?
      colwidth_++;
    }
    fmt_.append( "%" + left_arg + long_arg + width_arg + prec_arg + TypeChar_[type_] );
    colwidth_ += width_;
  }
}

bool TextFormat::IsDoubleType(FmtType typeIn) {
  return (typeIn == DOUBLE || typeIn == SCIENTIFIC || typeIn == GDOUBLE);
}

int TextFormat::SetFormatType(FmtType typeIn) {
  if (IsDoubleType(typeIn) && IsDoubleType(type_)) {
    type_ = typeIn;
    return 1;
  }
  return 0;
}

void TextFormat::SetCoordFormat(size_t maxFrames, double min, double step,
                                int default_width, int default_precision)
{
  int col_precision = default_precision;
  // Determine maximum coordinate.
  double maxCoord = (step * (double)maxFrames) + min;
  // Determine character width necessary to hold largest coordinate.
  int col_width = DigitWidth( (long int)maxCoord );
  // Check if the precision is enough to support the step size.
  if (step < 1.0) {
    int prec_exp_width = FloatWidth( step );
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
  if (type_ == INTEGER || type_ == STRING) // sanity check
    type_ = DOUBLE;
  width_ = col_width;
  precision_ = col_precision;
  align_ = RIGHT;
  SetFormatString();
}

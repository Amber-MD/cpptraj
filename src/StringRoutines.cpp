#include <cstdio>    // fprintf, fopen, fclose, sprintf
#include <cmath>     // log10
#include <sstream>   // istringstream, ostringstream
#include <locale>    // isspace
#include <stdexcept> // BadConversion
#include <vector>
#ifndef __PGI
#  include <glob.h>  // For tilde expansion
#endif

// tildeExpansion()
/** Use glob.h to perform tilde expansion on a filename, returning the 
  * expanded filename. The calling function is responsible for freeing
  * memory allocated with tildeExpansion.
  */
std::string tildeExpansion(const char *filenameIn) {
  if (filenameIn==0) {
    fprintf(stderr,"Error: tildeExpansion: null filename specified.\n");
    return std::string("");
  }
#ifdef __PGI
  // NOTE: It seems some PGI compilers do not function correctly when glob.h
  //       is included and large file flags are set. Just disable globbing
  //       for PGI and return a copy of filenameIn.
  return std::string(filenameIn);
#else
  glob_t globbuf;
  globbuf.gl_offs = 1;
  if ( glob(filenameIn, GLOB_TILDE, NULL, &globbuf)!=0 )
    return std::string("");
  //mprintf("DEBUG\tGLOB(0): [%s]\n", globbuf.gl_pathv[0]);
  std::string returnFilename( globbuf.gl_pathv[0] );
  globfree(&globbuf);
  return returnFilename;
#endif
}

// fileExists()
/** Return true if file can be opened "r".  */
bool fileExists(const char *filenameIn) {
  // Perform tilde expansion
  std::string fname = tildeExpansion(filenameIn);
  if (fname.empty()) return false;
  FILE *infile = fopen(fname.c_str(), "rb");
  if (infile==0) return false;
  fclose(infile);
  return true;
}

// NumberFilename()
/** Given a filename and a number, append number to filename, i.e.
  * filename.number.
  * The buffer should have enough space to handle the append.
  */
std::string NumberFilename(std::string const &fname, int number) {
  std::ostringstream oss;
  oss << fname << "." << number;
  return oss.str();
}

// DigitWidth()
/** \return the number of characters necessary to express the given digit. */
int DigitWidth(long int numberIn) {
  double numf;
  int minusSign = 0;

  if (numberIn == 0L) return 1;
  if (numberIn < 0L) {
    numf = (double)(-numberIn);
    minusSign = 1;
  } else
    numf = (double) numberIn;

  numf = log10( numf );
  ++numf;
  // The cast back to long int implicitly rounds down
  int numi = (int)numf;
  return (minusSign + numi);
}

// FloatWidth()
/** \return the number of characters necessary to express given float. */
int FloatWidth(double floatIn) {
  double float_exponent = fabs( log10( floatIn ) );
  ++float_exponent;
  return (int)float_exponent; // Cast to int implicitly rounds down
}

// ---------- STRING CONVERSION ROUTINES --------------------------------------- 
/*! \class: BadConversion
    \brief Runtime exception for catching bad conversions from the convertToX routines.
  */
class BadConversion : public std::runtime_error {
public:
  BadConversion(std::string const &s)
    : std::runtime_error(s)
    { }
};

// convertToInteger()
/// Convert the input string to an integer.
int convertToInteger(std::string const &s) {
  std::istringstream iss(s);
  long int i;
  if (!(iss >> i))
    throw BadConversion("convertToInteger(\"" + s + "\")");
  return (int)i;
}

// convertToDouble()
/// Convert the input string to a double.
double convertToDouble(std::string const &s) {
  std::istringstream iss(s);
  double d;
  if (!(iss >> d))
    throw BadConversion("convertToDouble(\"" + s + "\")");
  return d;
}

// RemoveTrailingWhitespace()
/// Remove any trailing whitespace from string.
void RemoveTrailingWhitespace(std::string &line) {
  std::locale loc;

  std::string::iterator p = line.end();
  --p;
  for (; p != line.begin(); p--)
    if (!isspace( *p, loc)) break;
  size_t lastSpace = (size_t)(p - line.begin()) + 1;
  //mprintf("lastSpace = %zu\n",lastSpace);
  if (lastSpace==1)
    line.clear();
  else
    line.resize( lastSpace );
}

std::string integerToString(int i) {
  std::ostringstream oss;
  oss << i;
  return oss.str();
}

std::string integerToString(int i, int width) {
  std::ostringstream oss;
  oss.fill('0');
  oss.width( width );
  oss << std::right << i;
  return oss.str();
}

std::string doubleToString(double d) {
  std::ostringstream oss;
  oss << d;
  return oss.str();
}

// ---------- STRING FORMAT ROUTINES -------------------------------------------
// SetDoubleFormatString()
/** Set up a printf-style format string for float/double of given width, 
  * precision, and alignment, e.g. '%8.3lf'.
  */
std::string SetDoubleFormatString(int width, int precision, int type)
{
  std::string format;
  std::string width_arg;
  std::string prec_arg;
  std::string type_arg; // Will be f, lf, or E.

  // Type: 1 = float, 2 = scientific (E), otherwise double
  switch (type) {
    case 1:  type_arg = "f"; break;
    case 2:  type_arg = "E"; break;
    default: type_arg = "lf"; break;
  }
  // Set width and/or precision if applicable.
  if (width > 0)
    width_arg = integerToString( width );
  if (precision > -1)
    prec_arg = "." + integerToString( precision );
  // Set format string.
  format.append( "%" + width_arg + prec_arg + type_arg );
  return format; 
}

// SetStringFormatString()
/** Set up a printf-style format string for string (char*) of given
  * width and alignment, e.g. '%20s'.
  */
std::string SetStringFormatString(int width, bool leftAlign)
{
  std::string format;
  std::string width_arg;
  // If not left-aligned, need leading space.
  if (!leftAlign) 
    format.assign("%");
  else
    format.assign("%-");
  // Set width if applicable
  if (width > 0)
    width_arg = integerToString( width );
  // Set format string.
  format.append( width_arg + "s" );
  return format;
}

// SetIntegerFormatString()
/** Set up a printf-style format string for integer of given width
  * and alignment, e.g. '%8i'.
  */
std::string SetIntegerFormatString(int width)
{
  std::string format;
  std::string width_arg;
  // Set width if applicable
  if (width > 0)
    width_arg = integerToString( width );
  // Set format string.
  format.append( "%" + width_arg + "i" );
  return format;
}

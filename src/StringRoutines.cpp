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
/** Return the number of characters necessary to express the given digit. */
int DigitWidth(int numberIn) {
  float numf;
  int numi;
  int minusSign = 0;

  if (numberIn==0) return 1;
  if (numberIn<0) {
    numf = (float) (-numberIn);
    minusSign = 1;
  } else
    numf = (float) numberIn;

  numf = log10( numf );
  ++numf;
  // The cast back to int implicitly rounds down
  numi = (int) numf;
  return (minusSign + numi);
}

// ---------- STRING FORMAT ROUTINES -------------------------------------------
// NOTE: In the following format routines, 2 char arrays are used
// since a 1 char null terminates the format string.

// SetDoubleFormatString()
/** Set up a printf-style format string for float/double of given width, 
  * precision, and alignment.
  */
void SetDoubleFormatString(std::string &formatString, int width, int precision,
                           int type, bool leftAlign)
{
  char leftSpace[2];
  char typestring[3];
  // If left-aligned, no leading space.
  if (leftAlign)
    leftSpace[0]='\0';
  else {
    leftSpace[0]=' ';
    leftSpace[1]='\0';
  }
  // Type: 1 = float, 2 = scientific (E), otherwise double
  if (type == 1) {
    typestring[0]='f';
    typestring[1]='\0';
  } else if (type == 2) {
    typestring[0] = 'l';
    typestring[1] = 'E';
    typestring[2] = '\0';
  } else {
    typestring[0] = 'l';
    typestring[1] = 'f';
    typestring[2] = '\0';
  }
  // # chars necessary to hold width arg
  int wWidth = DigitWidth( width );
  // # chars necessary to hold precision arg
  int pWidth = DigitWidth( precision );
  // String fmt: "%w.plf\0"
  char *format = new char[ pWidth + wWidth + 6 ];
  sprintf(format, "%s%%%i.%i%s", leftSpace, width, precision, typestring);
  formatString.assign( format );
  //mprintf("DEBUG: Double Format string: [%s]\n",format);
  delete[] format;
}

// SetStringFormatString()
/** Set up a printf-style format string for string (char*) of given
  * width and alignment.
  */
void SetStringFormatString(std::string &formatString, int width, bool leftAlign)
{
  char leftSpace[2];
  char alignChar[2];
  // If left-aligned, no leading space, set alignment char
  if (leftAlign) {
    leftSpace[0]='\0';
    alignChar[0]='-';
    alignChar[1]='\0';
  } else {
    leftSpace[0]=' ';
    leftSpace[1]='\0';
    alignChar[0]='\0';
  }
  // # chars necessary to hold width arg
  int wWidth = DigitWidth( width );
  // String fmt: " %-ws"
  char *format = new char[ wWidth + 5 ];
  sprintf(format, "%s%%%s%is", leftSpace, alignChar, width);
  formatString.assign( format );
  //mprintf("DEBUG: String Format string: [%s]\n",format);
  delete[] format;
}

// SetIntegerFormatString()
void SetIntegerFormatString(std::string &formatString, int width, bool leftAlign)
{
  char leftSpace[2];
  // If left-aligned, no leading space.
  if (leftAlign)
    leftSpace[0]='\0';
  else {
    leftSpace[0]=' ';
    leftSpace[1]='\0';
  }
  // # chars necessary to hold width arg
  int wWidth = DigitWidth( width );
  // String fmt: " %wi"
  char *format = new char[ wWidth + 4 ];
  sprintf(format, "%s%%%ii", leftSpace, width);
  formatString.assign( format );
  //mprinterr("DEBUG: Integer Format string: [%s]\n",format);
  delete[] format;
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

#include <cmath>     // log10
#include <cctype>    // isspace, isdigit
#include <ctime>     // for TimeString()
#include <sstream>   // istringstream, ostringstream
#include <stdexcept> // BadConversion
#include "StringRoutines.h"
#include "CpptrajStdio.h"
#ifdef _MSC_VER
# include <windows.h>
#else
# include <unistd.h>
#endif
#ifdef __APPLE__
# include <sys/sysctl.h>
# include <mach/mach_host.h>
#endif
#ifdef _WIN32
# include <time.h>
#endif

/** \param fname Input string.
  * \param number Input number.
  * \return string '<fname>.<number>'
  */
std::string AppendNumber(std::string const &fname, int number) {
  std::ostringstream oss;
  oss << fname << "." << number;
  return oss.str();
}

/** \return 1 if S1 matches S2. 0 otherwise.
  * \param S1 String containing wildcards.
  * \param S2 String to match.
  * The following wildcards are supported:
  *   '*': Any number of characters, including none.
  *   '?': A single character.
  */
int WildcardMatch(std::string const& S1, std::string const& S2) {
  std::string::const_iterator c1 = S1.begin();
  std::string::const_iterator c2 = S2.begin();
  while ( c1 != S1.end() || c2 != S2.end() ) {
    if (c1 == S1.end() && c2 != S2.end()) return 0;
    if (*c1 == '*') {
      ++c1;
      if (c1 == S1.end()) return 1;
      bool match = false;
      while (c2 != S2.end()) {
        if (*c2 == *c1) {
          match = true;
          break;
        }
        ++c2;
      }
      if (!match) return 0;
    } else if (c2 == S2.end()) {
      return 0;
    } else if (*c1 == '?') {
      ++c1;
      ++c2;
    } else {
      if (*c1 != *c2) return 0;
      ++c1;
      ++c2;
    }
  }
  return 1;
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
  iss >> i;
  if (iss.fail())
    throw BadConversion("convertToInteger(\"" + s + "\")");
  return (int)i;
}

// convertToDouble()
/// Convert the input string to a double.
double convertToDouble(std::string const &s) {
  std::istringstream iss(s);
  double d;
  iss >> d;
  if (iss.fail())
    throw BadConversion("convertToDouble(\"" + s + "\")");
  return d;
}

// RemoveTrailingWhitespace()
/// Remove any trailing whitespace from string.
void RemoveTrailingWhitespace(std::string &line) {
  if (line.empty()) return;
  int p = (int)line.size() - 1;
  while (p > -1 && (isspace(line[p]) || line[p]=='\n' || line[p]=='\r'))
    --p;
  line.resize(p + 1);
}

std::string NoTrailingWhitespace(std::string const& line) {
  std::string duplicate(line);
  RemoveTrailingWhitespace(duplicate);
  return duplicate;
}

void RemoveLeadingWhitespace(std::string& line) {
  size_t pos = 0;
  while (pos < line.size() && isspace(line[pos])) pos++;
  if (pos > 0)
    line = line.substr(pos, line.size() - pos);
}

std::string NoLeadingWhitespace(std::string const& line) {
  std::string duplicate(line);
  RemoveLeadingWhitespace(duplicate);
  return duplicate;
}

/// Remove all whitespace from string.
void RemoveAllWhitespace(std::string& line) {
  if (line.empty()) return;
  std::string tmp( line );
  line.clear();
  for (std::string::const_iterator it = tmp.begin(); it != tmp.end(); ++it) {
    if (isspace(*it) || *it == '\n' || *it == '\r') continue;
    line += *it;
  }
}

/// \return Given string with all whitespace removed.
std::string NoWhitespace(std::string const& line) {
  if (line.empty()) return std::string("");
  std::string out;
  for (std::string::const_iterator it = line.begin(); it != line.end(); ++it) {
    if (isspace(*it) || *it == '\n' || *it == '\r') continue;
    out += *it;
  }
  return out;
}

/// \return Given string with N leading characters removed.
std::string RemoveLeadingChars(std::string const& strIn, int NtoRemove)
{
  if (strIn.empty()) return std::string();
  std::string out = strIn.substr(NtoRemove, std::string::npos);
  return out;
}

// integerToString()
std::string integerToString(int i) {
  std::ostringstream oss;
  oss << i;
  return oss.str();
}

// integerToString()
std::string integerToString(int i, int width) {
  std::ostringstream oss;
  oss.fill('0');
  oss.width( width );
  oss << std::right << i;
  return oss.str();
}

// doubleToString()
std::string doubleToString(double d) {
  std::ostringstream oss;
  oss << d;
  return oss.str();
}

// validInteger()
bool validInteger(std::string const &argument) {
  if (argument.empty()) return false;
  std::string::const_iterator c;
  if (argument[0]=='-' || argument[0]=='+') {
    c = argument.begin()+1;
    if (c == argument.end()) return false;
  } else
    c = argument.begin();
  for (; c != argument.end(); ++c)
    if (!isdigit(*c)) return false;
  return true;
}

// validDouble()
bool validDouble(std::string const& argument) {
  if (argument.empty()) return false;
  // NOTE: We cannot just use an istringstream extraction operator >>
  //       because this will read *any* numbers that are there, so e.g.
  //       10cat is valid and gets converted to 10; the 'cat' is skipped.
  std::string::const_iterator ptr = argument.begin();
  bool isValid = false;
  if ( *ptr == '.' ||
       *ptr == '+' ||
       *ptr == '-' ||
       isdigit(*ptr) )
  {
    isValid = true;
    // Start of a number
    bool hasExponent = false;
    bool decimal_point = (*ptr == '.');
    if (!isdigit(*ptr)) ++ptr;
    while ( ptr != argument.end() )
    {
      // Check the next character
      if (*ptr == '.') {
        // Decimal point.
        if (hasExponent)
          // Do not allow decimal point in exponent
          return false;
        if (decimal_point)
          // Already have 1 decimal point. Do not allow two.
          return false;
        decimal_point = true;
      } else if (*ptr == 'E' || *ptr == 'e') {
        // Exponent.
        if (hasExponent)
          // Already have exponent. Do not allow two.
          return false;
        // Check start of exponent
        std::string::const_iterator ptr1 = ptr + 1;
        if (*ptr1 == '+' || *ptr1 == '-')
          // Advance past sign
          ++ptr;
        hasExponent = true;
      } else if ( !isdigit(*ptr) )
        return false;
      // Advance to next char
      ++ptr;
    }
  }
  return isValid;
/*
  std::istringstream iss(argument);
  double val;
  iss >> val;
  return !(iss.fail());
*/
}

// -----------------------------------------------------------------------------
// NOTE: I think this serves as a great example of how printf syntax is way
//       easier than iostream stuff (same printf command is only 3 lines). -DRR
std::string TimeString() {
  time_t rawtime;
  time( &rawtime );
  struct tm timeinfo;
# ifdef _WIN32
  localtime_s( &timeinfo, &rawtime );
# else
  localtime_r( &rawtime, &timeinfo );
# endif
  std::ostringstream oss;
  oss.fill('0');
  oss.width(2);
  oss << std::right << timeinfo.tm_mon+1;
  oss.put('/');
  oss.width(2);
  oss << std::right << timeinfo.tm_mday;
  oss.put('/');
  oss.width(2);
  oss << std::right << timeinfo.tm_year%100;
  oss.put(' ');
  oss.width(2);
  oss << std::right << timeinfo.tm_hour;
  oss.put(':');
  oss.width(2);
  oss << std::right << timeinfo.tm_min;
  oss.put(':');
  oss.width(2);
  oss << std::right << timeinfo.tm_sec;
  return oss.str();
}

// -----------------------------------------------------------------------------
std::string ByteString(unsigned long long sizeInBytes, ByteType bt) {
  static const char* BytePrefix[] = { " kB", " MB", " GB", " TB", " PB", " EB" };
  unsigned int idx = 0;    // Index into BytePrefix
  unsigned long long base; // Base of the size classes
  if (bt == BYTE_BINARY)
    base = 1024UL; // BINARY
  else
    base = 1000UL; // DECIMAL
  unsigned long long den = base;        // Number to divide input size by; matches idx
  unsigned long long cut = base * base; // Next size class
  while (sizeInBytes >= cut) {
    ++idx;
    den *= base;
    if (idx == 5) break; // NOTE: Must be max value of BytePrefix array.
    cut *= base;
  }
  double newSize = (double)sizeInBytes / (double)den;
  std::ostringstream oss;
  oss.setf( std::ios::fixed, std::ios::floatfield );
  oss.precision(3);
  oss << newSize;
  return oss.str() + std::string( BytePrefix[idx] );
}

// -----------------------------------------------------------------------------
#ifdef __APPLE__
static long long TotalGlobalMemory() {
    int mib[] = {CTL_HW, HW_MEMSIZE};
    int64_t size = 0;
    size_t len = sizeof(size);
    if (sysctl(mib, 2, &size, &len, NULL, 0) == 0)
        return (long long) size;
    return 0ll;
}
#endif

static long long AvailableMemory() {
#ifdef _MSC_VER
  MEMORYSTATUS status;
  GlobalMemoryStatus(&status);
  if (status.dwLength != sizeof(status))
    return -1;
  return (long long)status.dwAvailPhys;
#elif defined(__APPLE__)
  mach_msg_type_number_t count = HOST_VM_INFO_COUNT;
  vm_statistics_data_t vmstat;
  if (KERN_SUCCESS == host_statistics(mach_host_self(), HOST_VM_INFO,
                                      (host_info_t)&vmstat, &count)) {
    double total = vmstat.wire_count + vmstat.active_count +
                   vmstat.inactive_count + vmstat.free_count;
    double free = vmstat.free_count / total; // fraction
    return (long long)(TotalGlobalMemory() * free);
  }
  return -1;
#elif defined(_SC_AVPHYS_PAGES) && defined(_SC_PAGE_SIZE)
  long pages = sysconf(_SC_AVPHYS_PAGES);
  long page_size = sysconf(_SC_PAGE_SIZE);
  if (pages < 0L || page_size < 0L) return -1;
  return (long long)(pages * page_size);
#else
  return -1;
#endif
}

std::string AvailableMemoryStr() {
  long long avail_in_bytes = AvailableMemory();
  if (avail_in_bytes < 0)
    return std::string("");
  else
    return ByteString(avail_in_bytes, BYTE_DECIMAL);
}

/** \return A string containing the given array converted to range expression.
  */
std::string ArrayToRangeExpression(std::vector<int> const& arrayIn, int offsetIn) {
  if (arrayIn.empty())
    return std::string("");
  if (arrayIn.size() == 1)
    return integerToString( arrayIn.front() + offsetIn );

  std::string out("");
  // Commas will only be printed for groups after the first
  int commaGroup = 0;
  unsigned int idx = 0;
  while (idx < arrayIn.size()) {
    unsigned int kdx = idx + 1;
    int delta = arrayIn[kdx] - arrayIn[kdx-1];
    while (delta == 1) {
      kdx++;
      if (kdx == arrayIn.size()) break;
      delta = arrayIn[kdx] - arrayIn[kdx-1];
    }
    //mprintf("DEBUG: idx= %u kdx= %u array[i]= %i array[k-1]= %i delta= %i\n",
    //        idx, kdx, arrayIn[idx], arrayIn[kdx-1], delta);
    if (delta <= 0) {
      mprinterr("Internal Error: ArrayToRangeExpression() requires arrays in increasing order.\n");
      return std::string("");
    }
    if (commaGroup > 0)
      out.append(",");
    commaGroup++;
    unsigned int jdx = kdx - 1;
    if (idx == jdx)
      out.append( integerToString( arrayIn[idx] + offsetIn ) );
    else
      out.append( integerToString( arrayIn[idx] + offsetIn ) + "-" +
                  integerToString( arrayIn[jdx] + offsetIn ) );
    idx = kdx;
  }
 
  return out;
}

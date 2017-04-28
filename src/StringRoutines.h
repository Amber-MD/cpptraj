#ifndef INC_STRINGROUTINES_H
#define INC_STRINGROUTINES_H
#include <string>
/*! \file StringRoutines.h
    \brief Collection of useful string routines.

    Any commonly used routines that make use of the string class. This
    includes converting to/from numbers, string modification, setting
    printf-type format strings, etc.
 */
/// Append '.<number>' to string, useful for e.g. appending number to file name.
std::string AppendNumber(std::string const &, int);
/// \return 1 if first string (containing wildcards *, ?) matches second.
int WildcardMatch(std::string const&, std::string const&);
/// \return number of characters needed to represent given digit.
int DigitWidth(long int);
/// \return number of characters needed to represent given floating point.
int FloatWidth(double);
/// Remove any trailing whitespace from string.
void RemoveTrailingWhitespace(std::string &);
/// \return string stripped of trailing whitespace.
std::string NoTrailingWhitespace(std::string const&);
/// Remove all whitespace from a string
void RemoveAllWhitespace(std::string&);
/// \return string with all whitespace removed.
std::string NoWhitespace(std::string const&);
/// Convert string to integer.
int convertToInteger(std::string const &);
/// Convert string to double.
double convertToDouble(std::string const &);
/// Convert integer to string
std::string integerToString(int);
/// Convert integer to string, pad given width with zeros if necessary.
std::string integerToString(int,int);
/// Convert double to string
std::string doubleToString(double);
/// \return true if given string represents a valid integer.
bool validInteger(std::string const&);
/// \return true if given string represents a valid floating point number,
bool validDouble(std::string const&);
/// \return the current date/time with format 'mm/dd/yy  hh:mm:ss'
std::string TimeString();
/// BYTE_BINARY: Sizes are based on 1024^X. BYTE_DECIMAL: Sizes are based on 1000^X.
enum ByteType { BYTE_BINARY = 0, BYTE_DECIMAL };
/// \return Given size in bytes as a human-readable string.
std::string ByteString(unsigned long long, ByteType);
/// \return String containing size of available physical memory.
std::string AvailableMemoryStr();
#endif

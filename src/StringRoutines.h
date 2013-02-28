#ifndef INC_STRINGROUTINES_H
#define INC_STRINGROUTINES_H
#include <string>
#include <vector>
/*! \file StringRoutines.h
    \brief Collection of useful string routines.

    Any commonly used routines that make use of the string class. This
    includes converting to/from numbers, string modification, setting
    printf-type format strings, etc.
 */
std::string tildeExpansion(const char *);
bool fileExists(const char *);

std::string NumberFilename(std::string const &, int);
int DigitWidth(int);
void SetDoubleFormatString(std::string &, int, int, int, bool);
void SetStringFormatString(std::string &, int, bool);
void SetIntegerFormatString(std::string &, int, bool);

int convertToInteger(std::string const &);
double convertToDouble(std::string const &);
void RemoveTrailingWhitespace(std::string &);
std::string integerToString(int);
std::string integerToString(int,int);
std::string doubleToString(double);
#endif

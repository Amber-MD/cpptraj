#ifndef INC_CPPTRAJSTDIO_H
#define INC_CPPTRAJSTDIO_H
#include <string>
/*! \file CpptrajStdio.h
    \brief Interface between Cpptraj and Stdio.

    This is a useful abstraction that allows cpptraj to use several 
    often-used functions from the CSTDLIB stdio library without having
    to worry about details. For example, the mprintf function ensures
    that during parallel runs messages are only printed to the master
    thread, etc.
 */
#define OUTPUTFRAMESHIFT 1 ///< Used for output in DataFile and some TrajFiles
void mflush();
void mprintf(const char *, ...);
void mprinterr(const char *, ...);
void rprintf(const char *, ...);
void rprinterr(const char *, ...);
//void printerr(const char *, const char *, ...);
//void printwar(const char *, const char *, ...);
std::string tildeExpansion(const char *);
bool fileExists(const char *);
//void NumberFilename(char *, char *, int);
std::string NumberFilename(std::string const &, int);
int DigitWidth(int);
void SetDoubleFormatString(std::string &, int, int, int, bool);
void SetStringFormatString(std::string &, int, bool);
void SetIntegerFormatString(std::string &, int, bool);

int convertToInteger(std::string const &);
double convertToDouble(std::string const &);
void RemoveTrailingWhitespace(std::string &);
std::string integerToString(int);
std::string doubleToString(double);
#endif

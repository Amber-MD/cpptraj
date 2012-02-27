#ifndef INC_CPPTRAJSTDIO_H
#define INC_CPPTRAJSTDIO_H
/*! \file CpptrajStdio.h
    \brief Interface between Cpptraj and Stdio.

    This is a useful abstraction that allows cpptraj to use several 
    often-used functions from the CSTDLIB stdio library without having
    to worry about details. For example, the mprintf function ensures
    that during parallel runs messages are only printed to the master
    thread, etc.
 */
void mflush();
void mprintf(const char *, ...);
void mprinterr(const char *, ...);
void rprintf(const char *, ...);
void rprinterr(const char *, ...);
void printerr(const char *, const char *, ...);
void printwar(const char *, const char *, ...);
char *tildeExpansion(char *, int);
bool fileExists(char *);
void NumberFilename(char *, char *, int);
int DigitWidth(int);
#define OUTPUTFRAMESHIFT 1 ///< Used for output in DataFile and some TrajFiles
//char *SetDoubleFormatString(int, int);
//char *SetAlignedDoubleFormatString(int, int);
//char *SetStringFormatString(int);
//char *SetIntegerFormatString(int);
//char *SetXYZFormatString(int , int );
#endif

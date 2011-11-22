// CpptrajStdio.h
#ifndef INC_CPPTRAJSTDIO_H
#define INC_CPPTRAJSTDIO_H
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
char *SetDoubleFormatString(int, int);
char *SetStringFormatString(int);
char *SetIntegerFormatString(int);
char *SetXYZFormatString(int , int );
#endif

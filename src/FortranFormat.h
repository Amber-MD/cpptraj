// FortranFormat
// This file contains routines that pertain to reading and writing files with
// fortran FLAG and FORMAT keywords (currently only Amber Topology files).
#include "PtrajFile.h"
#include "Name.h"
// Enumerated type for Fortran Format
enum FortranFormat {
  UNKNOWN_FFORMAT=0, F10I8, F5E16_8, F20a4, F12I6, F3I8
};
#define NUMFORTRANFORMAT 6
// FFSIZE: Combined size of %FLAG and %FORMAT lines (81 * 2)
#define FFSIZE 162
int GetFortranBufferSize(FortranFormat,int,int);
void *getFlagFileValues(PtrajFile *,const char*,int,int);
char *DataToFortranBuffer(char*,const char*,FortranFormat,int*,double*,NAME*,int);

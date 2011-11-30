#include <cstring>
#include <cstdio>
#include <cstdarg>
#include "FileIO.h"

#ifdef MPI
#  include "MpiRoutines.h"
#endif

// FileIO::Printf()
/** Take the formatted string and write it to file using Write.
  */
void FileIO::Printf(const char *format, ...) {
  char buffer[1024];
  va_list args;
  va_start(args, format);
  vsprintf(buffer,format,args);
  this->Write(buffer, sizeof(char), strlen(buffer));
  va_end(args);
}

// FileIO::Rank_printf()
/** When MPI, printf only for the specified rank. If no MPI, behaves just
  * like above Printf.
  */
void FileIO::Rank_printf(int rank, const char *format, ...) {
  char buffer[1024];
  va_list args;
  va_start(args, format);
  vsprintf(buffer,format,args);
#ifdef MPI
    if (worldrank==rank)
#endif
      this->Write(buffer, sizeof(char), strlen(buffer));
  va_end(args);
}


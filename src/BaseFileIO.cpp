#include "BaseFileIO.h"
#include <cstring>

#ifdef MPI
#  include "PtrajMpi.h"
#endif

// BaseFileIO
/* BaseFileIO::Printf()
 * Take the formatted string and write it to file using Write.
 */
void BaseFileIO::Printf(const char *format, ...) {
  char buffer[1024];
  va_list args;
  va_start(args, format);
  vsprintf(buffer,format,args);
  this->Write(buffer, sizeof(char), strlen(buffer));
  va_end(args);
}

/* BaseFileIO::Rank_printf()
 * When MPI, printf only for the specified rank. If no MPI, behaves just
 * like above Printf.
 */
void BaseFileIO::Rank_printf(int rank, const char *format, ...) {
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


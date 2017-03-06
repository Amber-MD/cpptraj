#include <cstdio>
#include <cstdarg>
#ifdef MPI
#  include "Parallel.h"
#endif

enum IO_LEVEL_TYPE {
  IO_ALL = 0,    // Normal output
  IO_SILENT,     // Suppress STDOUT output
  IO_STAY_SILENT // Suppress All STDOUT output forever.
};
/// Controls mprintf/rprintf output.
static IO_LEVEL_TYPE world_io_level_ = IO_ALL;
/// Controls mprinterr/rprinterr output.
static bool suppressErrorMsg_ = false;
/// Where normal output should be written.
static FILE* STDOUT_ = stdout;

// mflush()
/** Call flush on STDOUT only if this is the master thread */
void mflush() {
# ifdef MPI
  if (!Parallel::World().Master()) return;
# endif
  fflush(STDOUT_);
}

/** Print message to STDOUT even if IO_SILENT */
void loudPrintf(const char* format, ...) {
# ifdef MPI
  if (!Parallel::World().Master()) return;
# endif
  va_list args;
  va_start(args,format);
  vfprintf(STDOUT_,format,args);
  va_end(args);
}

/** Print message to STDERR even if suppressErrorMsg_ */
void loudPrinterr(const char *format, ...) {
# ifdef MPI
  if (!Parallel::World().Master()) return;
# endif
  va_list args;
  va_start(args,format);
  vfprintf(stderr,format,args);
  va_end(args);
}
#ifdef PARALLEL_DEBUG_VERBOSE
// -----------------------------------------------------------------------------
/** Master prints message to STDOUT, others to mpidebugfile. */
void mprintf(const char*format, ...) {
  va_list args;
  va_start(args,format);
  if (Parallel::World().Master()) {
    vfprintf(STDOUT_, format, args);
    vfprintf(Parallel::mpidebugfile_, format, args);
  } else
    vfprintf(Parallel::mpidebugfile_, format, args);
  va_end(args);
}

/** Master prints message to STDERR, others to mpidebugfile. */
void mprinterr(const char *format, ...) {
  va_list args;
  va_start(args,format);
  if (Parallel::World().Master()) {
    vfprintf(stderr,format,args);
    vfprintf(Parallel::mpidebugfile_, format, args);
  } else
    vfprintf(Parallel::mpidebugfile_, format, args);
  va_end(args);
}
// -----------------------------------------------------------------------------
#else
/** Print message to STDOUT only if this is the master thread */
void mprintf(const char *format, ...) {
  if (world_io_level_ == IO_ALL) {
#   ifdef MPI
    if (!Parallel::World().Master()) return;
#   endif
    va_list args;
    va_start(args,format);
    vfprintf(STDOUT_,format,args);
    va_end(args);
  }
}

/** Print message to STDERR only if this is the master thread */
void mprinterr(const char *format, ...) {
  if (suppressErrorMsg_) return;
# ifdef MPI
  if (!Parallel::World().Master()) return;
# endif
  va_list args;
  va_start(args,format);
  vfprintf(stderr,format,args);
  va_end(args);
}
#endif
// rprintf()
/** Print message to STDOUT for this worldrank */
void rprintf(const char *format, ...) {
  if (world_io_level_ == IO_ALL) {
    va_list args;
    va_start(args, format);
#   ifdef MPI
    char buffer[1024];
    int nc = sprintf(buffer, "[%i]\t", Parallel::World().Rank());
    nc += vsprintf(buffer + nc, format, args);
    fwrite(buffer, 1, nc, STDOUT_);
#   else
    vfprintf(STDOUT_,format,args);
#   endif
    va_end(args);
  }
}

// rprinterr()
/** Print message to STDERR for this worldrank */
void rprinterr(const char *format, ...) {
  va_list args;
  if (suppressErrorMsg_) return;
  va_start(args,format);
# ifdef MPI
  char buffer[1024];
  int nc = sprintf(buffer, "[%i]\t", Parallel::World().Rank());
  nc += vsprintf(buffer + nc, format, args);
  fwrite(buffer, 1, nc, stderr);
# else
  vfprintf(stderr,format,args);
# endif
  va_end(args);
}

/** Change status of STDOUT output as long as SuppressAllOutput has not
  * been called.
  * \param silentIn if true, silence STDOUT output, otherwise enable.
  */
void SetWorldSilent(bool silentIn) {
  if (world_io_level_ != IO_STAY_SILENT) {
    if (silentIn)
      world_io_level_ = IO_SILENT;
    else
      world_io_level_ = IO_ALL;
  }
}

/** Suppress all STDOUT output for the entire run. */
void SuppressAllOutput() { world_io_level_ = IO_STAY_SILENT; }

void SuppressErrorMsg(bool suppressIn) { suppressErrorMsg_ = suppressIn; }

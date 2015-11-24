#include <cstdio>
#include <cstdarg>
#ifdef MPI
#  include "Parallel.h"
#endif

static bool worldsilent = false; // If true suppress all mprintf output.
static bool supressErrorMsg = false; // If true supress all mprinterr output.

// mflush()
/** Call flush on STDOUT only if this is the master thread */
void mflush() {
# ifdef MPI
  if (!Parallel::World().Master()) return;
# endif
  fflush(stdout);
}

/** Print message to STDOUT even if worldsilent */
void loudPrintf(const char* format, ...) {
# ifdef MPI
  if (!Parallel::World().Master()) return;
# endif
  va_list args;
  va_start(args,format);
  vfprintf(stdout,format,args);
  va_end(args);
}

/** Print message to STDERR even if supressErrorMsg */
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
    vfprintf(stdout, format, args);
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
  if (worldsilent) return;
# ifdef MPI
  if (!Parallel::World().Master()) return;
# endif
  va_list args;
  va_start(args,format);
  vfprintf(stdout,format,args);
  va_end(args);
}

/** Print message to STDERR only if this is the master thread */
void mprinterr(const char *format, ...) {
  if (supressErrorMsg) return;
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
  va_list args;
  if (worldsilent) return;
  va_start(args, format);
# ifdef MPI
  char buffer[1024];
  int nc = sprintf(buffer, "[%i]\t", Parallel::World().Rank());
  nc += vsprintf(buffer + nc, format, args);
  fwrite(buffer, 1, nc, stdout);
# else
  vfprintf(stdout,format,args);
# endif
  va_end(args);
}

// rprinterr()
/** Print message to STDERR for this worldrank */
void rprinterr(const char *format, ...) {
  va_list args;
  if (supressErrorMsg) return;
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

void SetWorldSilent(bool silentIn)   { worldsilent = silentIn;      }

void SupressErrorMsg(bool supressIn) { supressErrorMsg = supressIn; }

// printerr()
/** Print error message along with calling routine.  */
/*void printerr(const char *ROUTINE, const char *format, ...) {
  va_list args;

  va_start(args,format);
  fprintf(stdout,"Error: ");
  if (ROUTINE!=0)
    fprintf(stdout, "%s: ",ROUTINE);
  vfprintf(stdout,format,args);
  va_end(args);
  fprintf(stdout,"\n");
}*/

// printwar()
/** Print warning message along with calling routine.  */
/*void printwar(const char *ROUTINE, const char *format, ...) {
  va_list args;

  va_start(args,format);
  fprintf(stdout, "Warning: ");
  if (ROUTINE!=0)
    fprintf(stdout,"%s: ",ROUTINE);
  vfprintf(stdout,format,args);
  va_end(args);
  fprintf(stdout,"\n");
}*/

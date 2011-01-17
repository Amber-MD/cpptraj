/*
 * CpptrajStdio
 * Interface between Cpptraj and stdio.
 * May want to print messages in a parallel environment but dont want
 * to include all the parallel functionality.   
 * Also provide nice wrapping for warnings and error messages.
 * NOTE: Add ifdefs around worldrank stuff?
 */
#include <cstdio>
#include <cstdarg>
#ifdef MPI
#  include "PtrajMpi.h"
#endif

/*
 * mprintf()
 * Print message to STDOUT only if this is the master thread
 */
void mprintf(const char *format, ...) {
  va_list args;

#ifdef MPI
  if (worldrank!=0) return;
#endif
  va_start(args,format);
  vfprintf(stdout,format,args);
  va_end(args);
}

/*
 * mprinterr()
 * Print message to STDERR only if this is the master thread
 */
void mprinterr(const char *format, ...) {
  va_list args;

#ifdef MPI
  if (worldrank!=0) return; 
#endif
  va_start(args,format);
  vfprintf(stderr,format,args);
  va_end(args);
}

/*
 * rprintf()
 * Print message to STDOUT for this worldrank
 */
void rprintf(const char *format, ...) {
  va_list args;

#ifdef MPI
  fprintf(stdout,"[%i] ",worldrank);
#endif
  va_start(args,format);
  vfprintf(stdout,format,args);
  va_end(args);
  return;
}

/*
 * printerr()
 * Print error message along with calling routine.
 */
void printerr(const char *ROUTINE, const char *format, ...) {
  va_list args;

  va_start(args,format);
  fprintf(stdout,"Error: ");
  if (ROUTINE!=NULL)
    fprintf(stdout, "%s: ",ROUTINE);
  vfprintf(stdout,format,args);
  va_end(args);
  fprintf(stdout,"\n");
  return;
}

/*
 * printwar()
 * Print warning message along with calling routine.
 */
void printwar(const char *ROUTINE, const char *format, ...) {
  va_list args;

  va_start(args,format);
  fprintf(stdout, "Warning: ");
  if (ROUTINE!=NULL)
    fprintf(stdout,"%s: ",ROUTINE);
  vfprintf(stdout,format,args);
  va_end(args);
  fprintf(stdout,"\n");
  return;
}


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
#include <cstdlib> // tildeExpansion
#include <cstring> // tildeExpansion
#ifndef __PGI
#  include <glob.h> // For tilde expansion
#endif
#ifdef MPI
#  include "MpiRoutines.h"
#endif

/* mflush()
 * Call flush on STDOUT only if this is the master thread
 */
void mflush() {
#ifdef MPI
  if (worldrank!=0) return;
#endif
  fflush(stdout);
}

/* mprintf()
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

/* mprinterr()
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

/* rprintf()
 * Print message to STDOUT for this worldrank
 */
void rprintf(const char *format, ...) {
  va_list args;

#ifdef MPI
  fprintf(stdout,"[%i]\t",worldrank);
#endif
  va_start(args,format);
  vfprintf(stdout,format,args);
  va_end(args);
  return;
}

/* rprinterr()
 * Print message to STDERR for this worldrank
 */
void rprinterr(const char *format, ...) {
  va_list args;

#ifdef MPI
  fprintf(stderr,"[%i]\t",worldrank);
#endif
  va_start(args,format);
  vfprintf(stderr,format,args);
  va_end(args);
  return;
}

/* printerr()
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

/* printwar()
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

/* tildeExpansion()
 * Use glob.h to perform tilde expansion on a filename, returning the 
 * expanded filename. The calling function is responsible for freeing
 * memory allocated with tildeExpansion.
 */
char *tildeExpansion(char *filenameIn, int debug) {
  char *returnFilename;
#ifdef __PGI
  // NOTE: It seems some PGI compilers do not function correctly when glob.h
  //       is included and large file flags are set. Just disable globbing
  //       for PGI and return a copy of filenameIn.
  returnFilename = (char*) malloc( (strlen(filenameIn)+1) * sizeof(char));
  strcpy(returnFilename, filenameIn);
  return returnFilename;
#else
  glob_t globbuf;
  if (filenameIn==NULL) {
    mprinterr("Error: tildeExpansion: NULL filename specified.\n");
    return NULL;
  }
  globbuf.gl_offs = 1;
  if ( glob(filenameIn, GLOB_TILDE, NULL, &globbuf)!=0 ) return NULL;
  if (debug>1) mprintf("\tGLOB(0): [%s]\n",globbuf.gl_pathv[0]);
  returnFilename=(char*) malloc( (strlen(globbuf.gl_pathv[0])+1) * sizeof(char));
  strcpy(returnFilename, globbuf.gl_pathv[0]);
  globfree(&globbuf);
  return returnFilename;
#endif
} 

/* fileExists()
 * Return true if file can be opened "r".
 */
bool fileExists(char *filenameIn) {
  FILE *infile = NULL;
  char *fname;

  // Perform tilde expansion
  fname = tildeExpansion(filenameIn,0);
  if (fname==NULL) return false;
  infile=fopen(filenameIn,"rb");
  free(fname);
  if (infile==NULL) return false;
  fclose(infile);
  return true;
}

/* NumberFilename()
 * Given a filename and a number, append number to filename, i.e.
 * filename.number.
 * The buffer should have enough space to handle the append.
 */
void NumberFilename(char *buffer, char *filenameIn, int number) {
  sprintf(buffer,"%s.%i",filenameIn,number);
}


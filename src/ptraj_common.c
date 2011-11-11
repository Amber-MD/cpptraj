#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

// ----- Originally from utility.c -----
// Print functions: error and warning
// error()
void error(char *function, char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  fprintf(stderr, "\nERROR in %s: ", function);
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\n" );
  va_end(args);
  //exit(1);
}
// warning()
void warning(char *function, char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  fprintf(stderr, "\nWARNING in %s: ", function);
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\n" );
  va_end(args);
}

// Memory functions: safe_malloc, safe_realloc, and safe_free
// safe_malloc
void *safe_malloc(size_t size) {
  void *mem;
  if ( size ) {
    if ( (mem = (void *) malloc( size )) == NULL )
      error( "safe_malloc", "Error in alloc of %i bytes", size);
    mem = memset(mem, (char) 0, size);
  } else
    mem = NULL;
  return( mem );
}
// safe_realloc()
void *safe_realloc(void *mem, size_t cur_size, size_t increase) {
  void *temp;

  if ( cur_size == 0 )
    return ( safe_malloc( increase ) );
  else if ( (temp = (void *)
             realloc((void *) mem, (size_t) (cur_size + increase))) == NULL )
    error( "safe_realloc", "Couldn't allocate %i bytes more", increase);

  /* cast temp to avoid pointer arithmetic on void* which has unknown size */
  /* use char* as a byte level pointer per traditional malloc/realloc usage */
  memset( (char *) temp + cur_size, (char) 0, increase);

  mem = temp;
  return(mem);
}
// safe_free()
void safe_free(void *pointer) {
 if ( pointer != NULL ) free(pointer);
}

// FILE IO: safe_fopen and safe_fclose
// ----- Originally from io.c -----
FILE *safe_fopen(char *buffer, char *mode) {
  return( fopen(buffer, mode) );
} 
void safe_fclose(FILE *fileIn) {
  fclose( fileIn );
}


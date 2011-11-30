// FileIO_Mpi
#include <cstdio>
#include <cstdlib>
#include "FileIO_Mpi.h"

// CONSTRUCTOR
FileIO_Mpi::FileIO_Mpi() {
  pfile = (parallelType) malloc( sizeof(parallelType));
}

// DESTRUCTOR
FileIO_Mpi::~FileIO_Mpi() {
  free(pfile);
}

// FileIO_Mpi::Open()
int FileIO_Mpi::Open(const char *filename, const char *mode) {
  int err=0;

  switch( mode[0] ) {
    case 'r' : err=parallel_openFile_read(pfile, filename); break;
    case 'w' : err=parallel_open_file_write(pfile, filename); break;
    case 'a' : err=1; break; // NOTE: No MPI append for now
    default  : err=1; break;
  }

  return err;
}

// FileIO_Mpi::Close()
int FileIO_Mpi::Close() {
  parallel_closeFile(pfile);
  return 0;
}

// FileIO_Mpi::Read()
int FileIO_Mpi::Read(void *buffer, size_t size, size_t count) {
  //size_t numread;
  int numread;
  // Should never be able to call Read when fp is NULL.
  //if (fp==NULL) {
  //  fprintf(stdout,"Error: FileIO_Mpi::Read: Attempted to read NULL file pointer.\n");
  //  return 1;
  //}
  numread = parallel_fread(pfile, buffer, size * count);
  if (numread == -1) return -1;

  // NOTE: Check for errors here.
  return numread;
}

// FileIO_Mpi::Write()
int FileIO_Mpi::Write(void *buffer, size_t size, size_t count) {
  //size_t numwrite;
  // Should never be able to call Write when fp is NULL.
  //if (fp==NULL) {
  //  fprintf(stdout,"Error: FileIO_Mpi::Write: Attempted to write to NULL file pointer.\n");
  //  return 1;
  //}
  if ( parallel_fwrite(pfile, buffer, size * count) ) return 1;

  // NOTE: Check for errors here.
  return 0;
}

// FileIO_Mpi::Seek()
int FileIO_Mpi::Seek(off_t offset) {

  if ( parallel_fseek(pfile, offset, SEEK_SET) ) return 1;

  return 0;
}

// FileIO_Mpi::Rewind()
int FileIO_Mpi::Rewind() {
  if ( parallel_fseek(pfile, 0L, SEEK_SET) ) return 1;
  return 0;
}

/*long int FileIO_Mpi::Tell() {
  z_off_t zipOffset;

  zipOffset = gztell(fp);
  return (long int) zipOffset;
}*/

// FileIO_Mpi::Gets()
int FileIO_Mpi::Gets(char *str, int num) {

  if ( parallel_fgets(pfile,str,num) == NULL ) return 1;
  return 0;
}

// FileIO_Mpi::SetSize()
/** Set size of mpi file, required when splitting up writes.
  */
int FileIO_Mpi::SetSize(long int offset) {

  if ( parallel_setSize(pfile, offset) ) return 1;
  return 0;
}


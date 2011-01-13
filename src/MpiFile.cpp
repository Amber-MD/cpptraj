// MpiFile
#include <cstdio>
#include <cstdlib>
#include "MpiFile.h"

MpiFile::MpiFile() {
  pfile=(parallelType) malloc(sizeof(parallelType));
}

MpiFile::~MpiFile() {
  free(pfile);
}

int MpiFile::Open(const char *filename, const char *mode) {
  int err;

  switch( mode[0] ) {
    case 'r' : err=parallel_openFile_read(pfile, filename); break;
    case 'w' : err=parallel_open_file_write(pfile, filename); break;
    case 'a' : err=1; break; // NOTE: No MPI append for now
  }

  return 0;
}

int MpiFile::Close() {
  parallel_closeFile(pfile);
  return 0;
}

int MpiFile::Read(void *buffer, size_t size, size_t count) {
  //size_t numread;
  int numread;
  // Should never be able to call Read when fp is NULL.
  //if (fp==NULL) {
  //  fprintf(stdout,"Error: GzipFile::Read: Attempted to read NULL file pointer.\n");
  //  return 1;
  //}
  numread = parallel_fread(pfile, buffer, size * count);
  if (numread == -1) return -1;

  // NOTE: Check for errors here.
  return numread;
}

int MpiFile::Write(void *buffer, size_t size, size_t count) {
  //size_t numwrite;
  // Should never be able to call Write when fp is NULL.
  //if (fp==NULL) {
  //  fprintf(stdout,"Error: GzipFile::Write: Attempted to write to NULL file pointer.\n");
  //  return 1;
  //}
  if ( parallel_fwrite(pfile, buffer, size * count) ) return 1;

  // NOTE: Check for errors here.
  return 0;
}

int MpiFile::Seek(off_t offset) {

  if ( parallel_fseek(pfile, offset, SEEK_SET) ) return 1;

  return 0;
}

int MpiFile::Rewind() {
  if ( parallel_fseek(pfile, 0L, SEEK_SET) ) return 1;
  return 0;
}

/*long int MpiFile::Tell() {
  z_off_t zipOffset;

  zipOffset = gztell(fp);
  return (long int) zipOffset;
}*/

int MpiFile::Gets(char *str, int num) {

  if ( parallel_fgets(pfile,str,num) == NULL ) return 1;
  return 0;
}

/*
 * MpiFile::SetSize()
 * Set size of mpi file, required when splitting up writes.
 */
int MpiFile::SetSize(long int offset) {

  if ( parallel_setSize(pfile, offset) ) return 1;
  return 0;
}


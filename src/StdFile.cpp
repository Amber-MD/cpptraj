// StdFile: Standard C file operations
#include "StdFile.h" // BaseFileIO.h, cstdio

// CONSTRUCTOR
StdFile::StdFile() {
  fp = NULL;
  isStdout=false;
}

// DESTRUCTOR
StdFile::~StdFile() {
  if (fp!=NULL) this->Close();
}

/*
 * StdFile::Open()
 * Open file using standard C routines. If mode is WRITE and no
 * filename given default to stdout.
 */
int StdFile::Open(const char *filename, const char *mode) {
  if (filename==NULL) {
    if (mode[0]=='w') 
      fp=stdout;
    else
      return 1;
    isStdout=true;
  } else
    fp = fopen(filename, mode);
  if (fp==NULL) return 1;
  return 0;
}

/*
 * StdFile::Close()
 * Close stream if not stdout
 */
int StdFile::Close() {
  if (fp!=NULL && !isStdout) fclose(fp);
  fp=NULL;
  return 0;
}

/*
 * StdFile::Read()
 */
int StdFile::Read(void *buffer, size_t size, size_t count) {
  size_t numread;
  // Should never be able to call Read when fp is NULL.
  //if (fp==NULL) {
  //  fprintf(stdout,"Error: StdFile::Read: Attempted to read NULL file pointer.\n");
  //  return 1;
  //}
  numread = fread(buffer, size, count, fp);
  // NOTE: Check for errors here?
  if (ferror(fp)) {
    perror("Error during StdFile::Read");
    return -1;
  }
  if (feof(fp)) return -1;
  //if (numread!=(size*count)) return 1;
  return (int) numread;
}

/*
 * StdFile::Write()
 */
int StdFile::Write(void *buffer, size_t size, size_t count) {
  size_t numwrite;
  // Should never be able to call Write when fp is NULL.
  //if (fp==NULL) {
  //  fprintf(stdout,"Error: StdFile::Write: Attempted to write to NULL file pointer.\n");
  //  return 1;
  //}
  // DEBUG
  //char *temp;
  //temp=(char*) buffer;
  //printf("Calling standard write(%i): [%s]\n",size * count,temp);

  numwrite = fwrite(buffer, size, count, fp);
  // NOTE: Check for errors here.
  if (numwrite!=(size*count)) return 1;
  return 0;
}

/*
 * StdFile::Seek()
 * NOTE: Use fseeko for better compatibility with large files.
 */
int StdFile::Seek(off_t offset) {
  // DEBUG
  //printf("Calling standard seek(%i): %li\n",origin,offset);

  return fseeko(fp, offset, SEEK_SET);
}

/*
 * StdFile::Rewind()
 */
int StdFile::Rewind() {
  rewind(fp);
  return 0;
}

/*
 * StdFile::Tell()
 */
off_t StdFile::Tell() {
  return ftello(fp);
}

/*
 * StdFile::Gets()
 */
int StdFile::Gets(char *str, int num) {
  if ( fgets(str,num,fp) == NULL ) {
    //fprintf(stdout,"DEBUG: StdFile::Gets returned NULL (%s) %i\n",str,num);
    return 1;
  } else
    return 0;
}


// FileIO_Std: Standard C file operations
#include "FileIO_Std.h" // FileIO.h, cstdio

// CONSTRUCTOR
FileIO_Std::FileIO_Std() {
  fp = NULL;
  isStdout=false;
}

// DESTRUCTOR
FileIO_Std::~FileIO_Std() {
  if (fp!=NULL) this->Close();
}

// FileIO_Std::Open()
/** Open file using standard C routines. If mode is WRITE and no
  * filename given default to stdout.
  */
int FileIO_Std::Open(const char *filename, const char *mode) {
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

// FileIO_Std::Close()
/** Close stream if not stdout
  */
int FileIO_Std::Close() {
  if (fp!=NULL && !isStdout) fclose(fp);
  fp=NULL;
  return 0;
}

// FileIO_Std::Read()
int FileIO_Std::Read(void *buffer, size_t size, size_t count) {
  size_t numread;
  // Should never be able to call Read when fp is NULL.
  //if (fp==NULL) {
  //  fprintf(stdout,"Error: FileIO_Std::Read: Attempted to read NULL file pointer.\n");
  //  return 1;
  //}
  if (feof(fp)) return -1;
  numread = fread(buffer, size, count, fp);
  // NOTE: Check for errors here?
  if (ferror(fp)) {
    perror("Error during FileIO_Std::Read");
    return -1;
  }
  //if (numread!=(size*count)) return 1;
  return (int) numread;
}

// FileIO_Std::Write()
int FileIO_Std::Write(void *buffer, size_t size, size_t count) {
  size_t numwrite;
  // Should never be able to call Write when fp is NULL.
  //if (fp==NULL) {
  //  fprintf(stdout,"Error: FileIO_Std::Write: Attempted to write to NULL file pointer.\n");
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

// FileIO_Std::Seek()
// NOTE: Use fseeko for better compatibility with large files.
int FileIO_Std::Seek(off_t offset) {
  // DEBUG
  //printf("Calling standard seek(%i): %li\n",origin,offset);
#ifdef _MSC_VER
	return _fseeki64(fp,offset,SEEK_SET);
#else
  return fseeko(fp, offset, SEEK_SET);
#endif
}

// FileIO_Std::Rewind()
int FileIO_Std::Rewind() {
  rewind(fp);
  return 0;
}

// FileIO_Std::Tell()
off_t FileIO_Std::Tell() {
#ifdef _MSC_VER
	return _ftelli64(fp);
#else
  return ftello(fp);
#endif
}

// FileIO_Std::Gets()
int FileIO_Std::Gets(char *str, int num) {
  if ( fgets(str,num,fp) == NULL ) {
    //fprintf(stdout,"DEBUG: FileIO_Std::Gets returned NULL (%s) %i\n",str,num);
    return 1;
  } else
    return 0;
}


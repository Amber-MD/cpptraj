// FileIO_Std: Standard C file operations
#include <algorithm> // std::min, std::max
#include "FileIO_Std.h" // FileIO.h, cstdio

// CONSTRUCTOR
FileIO_Std::FileIO_Std() : fp_(NULL), isStream_(false) {}

// DESTRUCTOR
FileIO_Std::~FileIO_Std() { Close(); }

/** Open the specified stream. */
int FileIO_Std::OpenStream(StreamType type) {
  Close();
  switch (type) {
    case STDIN : fp_ = stdin; break;
    case STDOUT: fp_ = stdout; break;
    case STDERR: fp_ = stderr; break;
  }
  isStream_ = true;
  return 0;
}

// FileIO_Std::Open()
/** Open file using standard C routines. */
int FileIO_Std::Open(const char *filename, const char *mode) {
  if (filename == 0) return 1;
  Close();
  fp_ = fopen(filename, mode);
  if (fp_==NULL) return 1;
  isStream_ = false;
  return 0;
}

// FileIO_Std::Close()
/** Close file if not stream. */
int FileIO_Std::Close() {
  if (fp_!=NULL && !isStream_) fclose(fp_);
  fp_=NULL;
  isStream_ = false;
  return 0;
}

// FileIO_Std::Read()
int FileIO_Std::Read(void *buffer, size_t num_bytes) {
  size_t numread = fread(buffer, 1, num_bytes, fp_);
  if (ferror(fp_)) {
    perror("Error during FileIO_Std::Read");
    return -1;
  }
  return (int) numread;
}

// FileIO_Std::Write()
int FileIO_Std::Write(const void *buffer, size_t num_bytes) {
  size_t numwrite = fwrite(buffer, 1, num_bytes, fp_);
  // NOTE: Check for errors here.
  if (numwrite != num_bytes) return 1;
  return 0;
}

// FileIO_Std::Seek()
// NOTE: Use fseeko for better compatibility with large files.
int FileIO_Std::Seek(off_t offset) {
  // DEBUG
  //printf("Calling standard seek(%i): %li\n",origin,offset);
#ifdef _MSC_VER
	return _fseeki64(fp_,offset,SEEK_SET);
#else
  return fseeko(fp_, offset, SEEK_SET);
#endif
}

// FileIO_Std::Rewind()
int FileIO_Std::Rewind() {
  rewind(fp_);
  return 0;
}

// FileIO_Std::Tell()
off_t FileIO_Std::Tell() {
#ifdef _MSC_VER
	return _ftelli64(fp_);
#else
  return ftello(fp_);
#endif
}

// FileIO_Std::Gets()
int FileIO_Std::Gets(char *str, int num) {
  if ( fgets(str,num,fp_) == NULL ) {
    //fprintf(stdout,"DEBUG: FileIO_Std::Gets returned NULL (%s) %i\n",str,num);
    return 1;
  } else
    return 0;
}
